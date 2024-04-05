from src.rxn_ctr_mcs import *
from src.utils import load_json, save_json
from src.post_processing import *
from rdkit.Chem import AllChem
from collections import defaultdict
from src.post_processing import Enzyme, DatabaseEntry
from tqdm import tqdm

known_rxns = load_json("../data/mapping/known_rxns_swissprot_enzymes_plus_mcs_90_240310.json") # data/mapping/known_rxns_jni_w_enzyme_validation.json 
rule2krhash = defaultdict(list)
for k,v in known_rxns.items():

    # Convert enzymes and db entries to namedtuples
    enzymes = [Enzyme(*elt) for elt in v['enzymes']]
    v['enzymes'] = enzymes

    db_entries = [DatabaseEntry(*elt) for elt in v['db_entries']]
    v['db_entries'] = db_entries


kr_am_errors = {} # Track known rxn am errors
kekulize_issues = defaultdict(list)


am_ct = 0 # Number known rxns analyzed
kek_ct = 0 #
for n_kr, (krid, kr)  in tqdm(enumerate(known_rxns.items())):
    rxn_sma = kr['smarts']

    # Skip reactions that trigger RXNMapper atom mapping errors
    try:
        am_rxn_sma = atom_map(rxn_sma)
    except:
        am_ct += 1
        kr_am_errors[krid] = kr
        continue

    # Construct reaction object
    rxn = AllChem.ReactionFromSmarts(am_rxn_sma, useSmiles=True)
    rxn.Initialize()

    rc_atoms = rxn.GetReactingAtoms()

    # Construct rxn ctr mol objs
    for j, mol in enumerate(rxn.GetReactants()):
        try:
            rc = get_sub_mol(mol, rc_atoms[j])
        except:
            kekulize_issues[krid].append((Chem.MolToSmiles(mol), rc_atoms[j]))
            kek_ct += 1
            continue

save_json(kr_am_errors, "../artifacts/analysis_issues/atom_mapping_issues_known_rxns_swissprot_plus_mcs_90.json")
save_json(kekulize_issues, "../artifacts/analysis_issues/get_reaction_center_mol_rdkit_issues_known_rxns_swissprot_plus_mcs_90.json")
print(f"{len(kr_am_errors)} reactions couldn't be atom mapped, {len(kekulize_issues)} reaction centers had issues getting sub mol")