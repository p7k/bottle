from src.utils import load_json
from collections import defaultdict
from src.post_processing import Enzyme, DatabaseEntry

known_rxns = load_json("../data/mapping/known_rxns_swissprot_enzymes_plus_mcs_90_240310.json") # data/mapping/known_rxns_jni_w_enzyme_validation.json 
rule2krhash = defaultdict(list)
for k,v in known_rxns.items():

    # Convert enzymes and db entries to namedtuples
    enzymes = [Enzyme(*elt) for elt in v['enzymes']]
    v['enzymes'] = enzymes

    db_entries = [DatabaseEntry(*elt) for elt in v['db_entries']]
    v['db_entries'] = db_entries