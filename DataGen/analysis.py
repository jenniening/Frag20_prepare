import os
import rdkit
from rdkit import Chem

class getInfor:
    """
    Get SMILES and InChI:
    infile1: the initial SMILES
    infile2: the sdf file for MMFF optimized structure
    infile3: the sdf file for QM optimized structure
    infor will include all 6 informations
    """
    def __init__(self,infile1, infile2, infile3, isomericSmiles):
        self._infile = [infile1, infile2, infile3]
        self._ftype = ["smiles", "sdf", "sdf"]
        self._isomericSmiles = isomericSmiles
        self._infor = []
        self._inforname = ["SMILES", "initial_InChI", "initial_SMILES", "MMFF_InChI", "MMFF_SMILES", "QM_InChI", "QM_SMILES"]

    @property
    def infile(self):
        return self._infile
    
    @property
    def ftype(self):
        return self._ftype
    
    @property
    def isomericSmiles(self):
        return self._isomericSmiles
    
    @property
    def inforname(self):
        return self._inforname

    @property
    def infor(self):
        if self._infor == []:
            self.getAll3()
        return self._infor


    def getAll3(self):
        """
        Get information for all three files, if this process failed, the output will be None, None
        """
        self._infor.append(self._infile[0])
        for idx, infile in enumerate(self._infile):
            try:
                self._infor.extend(self.getInfor(infile, self._ftype[idx], True, self._isomericSmiles))
            except:
                self._infor.extend(["None","None"])
    
    def getInfor(self, infile, file_type, removeHs, isomericSmiles):
        """
        Generate SMILES and InChI strings
        """
        if file_type == "sdf":
            mol = Chem.SDMolSupplier(infile, removeHs=removeHs)[0]
        else:
            mol = Chem.MolFromSmiles(infile)
        Inchi = Chem.MolToInchi(mol)
        Smiles = Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles)
        return [Inchi, Smiles]


def convert_tosdf(datadir, index):
    """
    Convert G09 log to SDF file
    Input file: xxxx.opt.log
    Output file: xxxx.opt.sdf
    """
    fn = os.path.join(datadir, str(index)+ '.opt.log')
    outfn = os.path.join(datadir, str(index)+ '.opt.sdf')
    cmd = 'obabel -ig03 ' + fn + ' -osdf -O ' + outfn
    os.system(cmd)

if __name__ == "__main__":
    out = open("Infor_isomeric.csv", "w")
    index_file = "/beegfs/jl7003/ccdc/script/Index_ccdc_20_overlapwithmol20_removeempty.csv"
    error_file = "/beegfs/jl7003/ccdc/script/Error_20_iteration_1.csv"
    error_list = [line.rstrip() for line in open(error_file)]
    index_list = [line.split(",")[0] for line in open(index_file) if line.split(",")[0] not in error_list][1:]
    smiles = [line.split(",")[4].rstrip() for line in open(index_file) if line.split(",")[0] not in error_list][1:]
    datadir = "/beegfs/jl7003/ccdc/cry_min"
    isomericSmiles = True
    for idx, i in enumerate(index_list):
        ### first, convert log file into sdf file ###
        convert_tosdf(datadir, i)
        infile1 = smiles[idx]
        infile2 = os.path.join(datadir, str(i) + ".sdf")
        infile3 = os.path.join(datadir, str(i) + ".opt.sdf")
        infile = getInfor(infile1, infile2, infile3, isomericSmiles)
        out.write(str(i) + "," + infile_1 + "," + ",".join(infile.infor) + "\n")
    out.close()

    
