# Define classes for pathway and reaction entries
from src.utils import sort_x_by_y
from src.pathway_utils import get_stoich_pk
import PIL
from rdkit.Chem import Draw
from rdkit.Chem import AllChem


class pathway:
    def __init__(self, rhashes, starter_hash=None, target_hash=None, prc_mcs=None, dG=None):
        self.starter = starter_hash
        self.target = target_hash
        self.rhashes = rhashes # Hash ids for the path's reactions, in order
        self.prc_mcs = prc_mcs # Average (over known rxns) peri-rxn-ctr MCS score for each predicted rxn (tuple)
        self.dG = dG # Placeholder for thermo

    def min_mcs(self):
        if self.prc_mcs is None:
            return None
        else:
            return min(self.prc_mcs)
        
    def max_mcs(self):
        if self.prc_mcs is None:
            return None
        else:
            return max(self.prc_mcs)
        
    def mean_mcs(self):
        if self.prc_mcs is None:
            return None
        else:
            return sum(self.prc_mcs) / len(self.prc_mcs)

    def compute_mean_prc_mcs(self, pred_rxns):
        '''
        Assumes known reactions sorted by 
        substrate_averaged_prc_mcs and could
        break. Need to ensure entries always sorted
        in expansion object.
        '''
        self.prc_mcs = []
        for rh in self.rhashes:
            krs = pred_rxns[rh].known_rxns
            kr_mean_mcs = 0
            for i, elt in enumerate(krs):
                if elt[0] is None:
                    break
                elif i > 0:
                    kr_mean_mcs = (kr_mean_mcs * i + sum(elt[0]) / len(elt[0])) / (i + 1) # Rolling ave
                else:
                    kr_mean_mcs = sum(elt[0]) / len(elt[0])

            self.prc_mcs.append(kr_mean_mcs)

class reaction:
    def __init__(self, rid, smarts, rules=[], known_rxns=[]):
        self.rid = rid
        self.smarts = smarts
        self.rules = rules
        self.known_rxns = known_rxns

    def sort_known_rxns(self):
        '''
        Sort by mean for now. May include input to sort by min
        '''
        krs_w_mcs = [elt for elt in self.known_rxns if elt[0] is not None]
        krs_wo_mcs = [elt for elt in self.known_rxns if elt[0] is None]

        if krs_w_mcs:
            mcses = list(zip(*krs_w_mcs))[0]
            mean_mcses = list(map(lambda x: sum(x) / len(x), mcses))
            krs_w_mcs, _ = sort_x_by_y(krs_w_mcs, mean_mcses, reverse=True)
            self.known_rxns = list(krs_w_mcs) + krs_wo_mcs

# class expansion:
#     def __init__(self, paths=[], pred_rxns={}):
#         self.paths = paths
#         self.pred_rxns = pred_rxns

def rxn_hash_2_rxn_sma(rhash, pk):
    '''
    Make reaction smarts string for
    reaction indexed by rhash in a pk
    object
    '''
    rxn_stoich = get_stoich_pk(rhash, pk)
    products = ".".join([".".join([smi]*stoich) for smi, stoich in rxn_stoich.items() if stoich >= 0])
    reactants = ".".join([".".join([smi]*abs(stoich)) for smi, stoich in rxn_stoich.items() if stoich <= 0])
    rxn_sma = ">>".join([reactants, products])
    return rxn_sma

# Pathway drawing functions
def draw_rxn(rxn_sma):
    return Draw.ReactionToImage(
        AllChem.ReactionFromSmarts(rxn_sma, useSmiles=True),
        subImgSize=(200, 200), useSVG=False, drawOptions=None, returnPNG=False
    )

def get_concat_h(im1, im2):
    dst = PIL.Image.new('RGB', (im1.width + im2.width, max(im1.height, im2.height)))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    return dst

def get_concat_v(im1, im2):
    dst = PIL.Image.new('RGB', (max(im1.width, im2.width), im1.height + im2.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst

def draw_pathway(pred_known_pairs):
    for i, elt in enumerate(pred_known_pairs):
        if i == 0:
            img = get_concat_h(*elt)
        else:
            img = get_concat_v(img, get_concat_h(*elt))

    return img