import pandas as pd
from minedatabase.pickaxe import Pickaxe
from collections import defaultdict
from argparse import ArgumentParser
from multiprocessing import set_start_method

from src.config import filepaths
from src.post_processing import Path, PredictedReaction, KnownReaction, Enzyme, DatabaseEntry, get_path_id
from src.rcmcs import extract_operator_patts, calc_lhs_rcmcs
from src.operator_mapping import expand_paired_cofactors, expand_unpaired_cofactors, standardize_template_map
from src.pickaxe_processing import find_paths, prune_pickaxe
from src.utils import load_json, save_json
from src.chem_draw import draw_rxn_svg

from src.thermo.batch_add_eq_compounds import add_compounds_to_eQ
from equilibrator_assets.local_compound_cache import LocalCompoundCache
from equilibrator_cache.compound_cache import CompoundCache
import sqlalchemy
from src.thermo.pickaxe_thermodynamics import PickaxeThermodynamics

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("filename", help='Expansion filename not including the extension .pk', type=str)
    parser.add_argument("generations", help="Number of generations run in this expansion", type=int)
    parser.add_argument("--do_thermo", action="store_true", help="Does thermo calculations if provided")
    args = parser.parse_args()

    if args.do_thermo:
        set_start_method("spawn")
    
    # Set params
    k_tautomers = 10 # How many top scoring tautomers to generate for operator mapping
    pre_standardized = False # Predicted reactions assumed pre-standardized
    
    # CMD params
    generations = args.generations
    pk_path = filepaths['raw_expansions'] / f"{args.filename}.pk"

    # Load stored paths
    path_filepath = filepaths['processed_expansions'] / 'found_paths.json'
    predicted_reactions_filepath = filepaths['processed_expansions'] / "predicted_reactions.json"
    known_reactions_filepath = filepaths['processed_expansions'] / "known_reactions.json"
    load_processed = lambda path : load_json(path) if path.exists() else {}
    stored_paths = load_processed(path_filepath)
    stored_predicted_reactions = load_processed(predicted_reactions_filepath)
    stored_known_reactions = load_processed(known_reactions_filepath)

    # Read in rules
    rules_dir = filepaths['rules']
    rules_fns = ["minimal1224_all_uniprot.tsv", "JN3604IMT_rules.tsv"]
    read_pd = lambda fn : pd.read_csv(f"{rules_dir}/{fn}", sep='\t').set_index("Name").drop(columns='Comments')
    min_rules, imt_rules = [read_pd(fn) for fn in rules_fns]

    # Read in known reactions
    # NOTE: Minified IMT operator must match MIN operator of a known reaction to be counted
    known_reaction_bank = load_json(filepaths['data'] / "sprhea/sprhea_240310_v3_mapped_no_subunits.json")
    imt2krs = defaultdict(list)
    for k, v in known_reaction_bank.items():
        if v['imt_rules'] and v['min_rule']:
            for imt in v['imt_rules']:
                if imt.split('_')[0] == v['min_rule']:
                    imt2krs[imt].append(k)

    imt2ct = {k : len(v) for k, v in imt2krs.items()}

    # Read in cofactor lookup tables
    paired_ref = pd.read_csv(filepaths['cofactors'] / 'paired_cofactors_reference.tsv', sep='\t')
    unpaired_ref = pd.read_csv(filepaths['cofactors'] / 'unpaired_cofactors_reference.tsv', sep='\t')
    smi2paired_cof = expand_paired_cofactors(paired_ref, k=k_tautomers)
    smi2unpaired_cof = expand_unpaired_cofactors(unpaired_ref, k=k_tautomers)

    # Load raw expansion object
    pk = Pickaxe()
    pk.load_pickled_pickaxe(pk_path)

    print("Finding paths")
    paths, starters, targets = find_paths(pk, generations)

    pk = prune_pickaxe(pk, paths)
    print(f"Pruned pk object to {len(pk.compounds)} compounds and {len(pk.reactions)} reactions")

    if args.do_thermo:
        print("Adding compounds to equilibrator")
        add_compounds_to_eQ(pk)

    # Create new PredictedReaction objects where don't already have
    print("Creating new predicted reactions")
    new_predicted_reactions = {}
    for sid, tid in paths.keys():
        for path in paths[(sid, tid)]:        
            for rid in path:
                if rid not in stored_predicted_reactions:
                    pr = PredictedReaction.from_pickaxe(pk, rid)
                    new_predicted_reactions[rid] = pr

    # Add KnownReactions to new PredictedReactions
    print("Adding known analogues")
    new_known_reactions = {}
    bad_ops = []
    for id, pr in new_predicted_reactions.items():
        analogues = {}
        imt_that_mapped_krs = [elt for elt in pr.operators if elt in imt2ct] # Filter out those that don't map any known reactions
        srt_imt = sorted(imt_that_mapped_krs, key= lambda x : imt2ct[x], reverse=True) # If multiple imt operators, start w/ most common
        for imt in srt_imt:
            min = imt.split('_')[0] # Minify imt operator to get reaction center by protection-guess-and-check
                
            did_map, aligned_smarts, reaction_center = standardize_template_map(
                rxn=pr.smarts,
                rule_row=min_rules.loc[min],
                smi2paired_cof=smi2paired_cof,
                smi2unpaired_cof=smi2unpaired_cof,
                return_rc=True,
                pre_standardized=pre_standardized,
                quiet=True,
            )

            if did_map:
                pr.smarts = aligned_smarts
                pr.reaction_center = reaction_center
                lhs_patts = extract_operator_patts(min_rules.loc[min, 'SMARTS'], side=0)
                pr_rcts = pr.smarts.split(">>")[0].split('.')
                pr_rcts_rc = [pr_rcts, pr.reaction_center]

                for krid in imt2krs[imt]: # Assign analogue to pr only on imt operator level
                    if krid in stored_known_reactions: # Load from stored known reactions
                        kr = KnownReaction.from_dict(stored_known_reactions[krid])
                    else: # Create new known reaction from bank
                        bank_kr = known_reaction_bank[krid]
                        
                        # Combine all known reaction operators
                        combined_ops = []
                        if bank_kr['min_rule']:
                            combined_ops.append(bank_kr['min_rule'])
                        if bank_kr['imt_rules']:
                            combined_ops += bank_kr['imt_rules']
                        
                        # Create known reaction object
                        kr = KnownReaction(
                            id=krid,
                            smarts=bank_kr['smarts'],
                            operators=combined_ops,
                            enzymes=[Enzyme.from_dict(e) for e in bank_kr['enzymes']],
                            db_entries=[DatabaseEntry.from_dict({'name': 'rhea', 'id': rhea}) for rhea in bank_kr['rhea_ids']],
                            reaction_center=bank_kr['reaction_center'],
                        )
                        new_known_reactions[krid] = kr # Store in dict of new krs
                        
                    # RCMCS
                    kr_rcts_rc = [
                        kr.smarts.split('>>')[0].split('.'), # Reactants
                        kr.reaction_center, # Reaction center
                    ]
                    rcmcs = calc_lhs_rcmcs(pr_rcts_rc, kr_rcts_rc, patts=lhs_patts, norm='max')
                    pr.rcmcs[krid] = rcmcs

                    analogues[krid] = kr # Append predicted reaction analogues
            else:
                print(f"Minified operator failed to recapitulate reaction {imt} {id}")
                bad_ops.append((imt, id))

            pr.analogues = analogues # Add analogues to predicted reaction

    if args.do_thermo:
        # Connect to compound cache
        with open(filepaths['artifacts'] / "eq_uris.uri", "r") as f:
            URI_EQ = f.read().strip("\n")
        
        lcp = LocalCompoundCache()
        lcp.ccache = CompoundCache(sqlalchemy.create_engine(URI_EQ))
    
        # Create pk thermo and eQ objects
        print(f"Getting Thermo Values for {len(pk.compounds)} compounds and {len(pk.reactions)} reactions")
        PT = PickaxeThermodynamics(lc=lcp)
        PT.generate_eQ_compound_dict_from_pickaxe(pk=pk)
        PT.generate_eQ_reaction_dict_from_pickaxe(pk=pk)
    
    # Create Path objects
    print("Adding new paths (and calculating mdf)")
    new_paths = {}
    for sid, tid in paths.keys():
        for path in paths[(sid, tid)]:
            pid = get_path_id(path)

            prs = []
            for rid in path:
                if rid in new_predicted_reactions:
                    prs.append(new_predicted_reactions[rid])
                else:
                    prs.append(PredictedReaction.from_dict(stored_predicted_reactions[rid], stored_known_reactions))

            # If new path or missing thermo, and want to do thermo now, do it
            if (pid not in stored_paths or stored_paths[pid]['mdf'] is None) and args.do_thermo:
                mdf_res = PT.calculate_pathway_mdf(reaction_id_list=[pr.id for pr in prs])
                mdf = mdf_res.mdf_value
                dG_opt = {k : v.magnitude for k,v in mdf_res.reaction_energies.items()}
                dG_err = {k : v.magnitude for k,v in mdf_res.uncertainties.items()}
                
                if mdf is None:
                    print(f"Failed mdf for path: {pid}")
            
            # Otherwise these are the defaults
            else:
                mdf = None
                dG_opt = {}
                dG_err = {}

            # If it was new, create a new path
            if pid not in stored_paths:
                # Add new path
                new_paths[pid] = Path(
                    id=pid,
                    starter=starters[sid],
                    target=targets[tid],
                    reactions=prs,
                    mdf=mdf,
                    dG_opt=dG_opt,
                    dG_err=dG_err,
                    sid=sid,
                    tid=tid,
                )
            
            # If it was missing thermo and old, update stored paths
            elif stored_paths[pid]['mdf'] is None and args.do_thermo:
                stored_paths[pid]['mdf'] = mdf
                stored_paths[pid]['dG_opt'] = dG_opt
                stored_paths[pid]['dG_err'] = dG_err

    # Generate rxn svgs
    for prid, pr in new_predicted_reactions.items():
        pr.image = draw_rxn_svg(pr.smarts, pr.id)

    for krid, kr in new_known_reactions.items():
        kr.image = draw_rxn_svg(kr.smarts, kr.id)

    # Add new to old
    new = [new_known_reactions, new_predicted_reactions, new_paths]
    old = [stored_known_reactions, stored_predicted_reactions, stored_paths]

    for n, o in zip(new, old):
        for id, entry in n.items():
            o[id] = entry.to_dict()

    # Save
    save_json(stored_paths, path_filepath)
    save_json(stored_predicted_reactions, predicted_reactions_filepath)
    save_json(stored_known_reactions, known_reactions_filepath)