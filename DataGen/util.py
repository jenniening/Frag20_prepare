import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def get_RMSD(prob, ref):
    ''' Get the RMSD between local minimum and global minimum 
        prob:conformations
        ref: crystal or reference structure
    '''
    
    RMSD = None
    heavyatomidx = []
    for a in ref.GetAtoms():
            if a.GetAtomicNum() != 1:
                heavyatomidx.append(a.GetIdx())
    try:
        RMSD = Chem.rdMolAlign.GetBestRMS(prob,ref)
    except:
        RMSD = Chem.rdMolAlign.AlignMol(prob,ref,-1,-1,atomMap = [(k, k) for k in heavyatomidx])
    return RMSD

