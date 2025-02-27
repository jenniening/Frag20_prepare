import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import os

def get_RMSD(prob, ref):
    ''' Get the RMSD between local minimum and global minimum 

    Args:
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


def convert_logtosdf(datadir, index):
    """
    Convert G09 log to SDF file
    Input file: xxxx.opt.log
    Output file: xxxx.opt.sdf
    """
    fn = os.path.join(datadir, str(index)+ '.opt.log')
    outfn = os.path.join(datadir, str(index)+ '.opt.sdf')
    cmd = 'obabel -ig03 ' + fn + ' -osdf -O ' + outfn
    os.system(cmd)
