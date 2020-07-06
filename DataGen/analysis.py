import os
import rdkit
from rdkit import Chem
import rdkit.Chem.Descriptors as Descriptors
import pandas as pd
from DataGen.util import convert_logtosdf

class getInfor:
    """Get SMILES and InChI

    Args:
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
            self._getAll3()
        return self._infor


    def _getAll3(self):
        """
        Get information for all three files, if this process failed, the output will be None, None
        """
        self._infor.append(self._infile[0])
        for idx, infile in enumerate(self._infile):
            try:
                self._infor.extend(self._getInfor(infile, self._ftype[idx], True, self._isomericSmiles))
            except:
                self._infor.extend(["None","None"])
    
    def _getInfor(self, infile, file_type, removeHs, isomericSmiles):
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

class check:
    """
    Analysis data after QM calculations
    """
    def __init__(self, datadir, index_list, outdir, confs=False):
        """

        Args:
            datadir: directory for all data files including MMFF optimized sdf file and QM calculated log file
            outdir: directory for all output data index files 
            index_list: index_list for all successful calculated files. Exp. 1.opt.log, 2.opt.log --> index_list[1,2]
            confs: if calculated files are from different conformations, consf=True, else, False. Defaults to False
        """
        self.datadir = os.path.abspath(datadir)
        self.outdir = os.path.abspath(outdir)
        self.olddir = os.getcwd()
        self.index_list = index_list
        self.confs = confs
        if self.confs:
            self.suffix = ".sdf"
        else:
            self.suffix = "_min.sdf"
        self.smiles_list, self.index_init_list = self.__get_relation__()

    def __get_relation__(self):
        """
        Get initial index list and smiles list for each calculated file
        """
        os.chdir(self.datadir)
        smiles_list = []
        init_index_list = []
        for i in self.index_list:
            with open(str(i) + self.suffix) as infile:
                lines = infile.readlines()
                index_init = lines[-9].rstrip()
                smiles_init = lines[-15].rstrip()
                smiles_list.append(smiles_init)
                init_index_list.append(index_init)
        os.chdir(self.olddir)
        return smiles_list, init_index_list

    def update_list(self,infile):
        """Update current index_list and initial_index_list based on provided file

        Args:
            infile: data index file
        """
        data = pd.read_csv(os.path.join(self.outdir,infile))
        self.index_list = list(data["index"])
        if "initial_index" in data.columns:
            self.initial_index_list = list(data["initial_index"])

    def built_initialdata(self):
        """
        Build initial information data
        """
        os.chdir(self.datadir)
        out = open(os.path.join(self.outdir,"initial_dataset.csv"), "w")
        out.write("index initial_index SMILES initial_InChI initial_SMILES MMFF_InChI MMFF_SMILES QM_InChI QM_SMILES\n")
        for idx, i in enumerate(self.index_list):
            infile1 = self.smiles_list[idx]
            infile2 = str(i) + self.suffix
            infile3 = str(i) + ".opt.sdf"
            if not os.path.exists(infile3):
                convert_logtosdf(self.datadir, i)
            infile = getInfor(infile1, infile2, infile3, isomericSmiles=False)
            out.write(str(i) + " " + self.index_init_list[idx] +  " " + " ".join(infile.infor) + "\n")
        out.close()
        os.chdir(self.olddir)

    def check_consistency(self, rule):
        """
        Check consistency of initial data, MMFF optimized data, and QM optimized data. 
        Remove structures which are not consistent based on rule
        
        Args:
            rule: "strict" or "loose"
        """
        data = pd.read_csv(os.path.join(self.outdir, "initial_dataset.csv"), sep=" ")
        if rule == "strict":
            data = data[(data["initial_SMILES"] == data["MMFF_SMILES"]) & (data["MMFF_SMILES"] == data["QM_SMILES"])]
            data.to_csv(os.path.join(self.outdir, "data_consistent_strict.csv"), index = False)
        elif rule == "loose":
            data_initial = data[(data["initial_SMILES"] == data["MMFF_SMILES"]) | (data["initial_InChI"] == data["MMFF_InChI"])]
            data = data_initial[(data_initial["MMFF_SMILES"] == data_initial["QM_SMILES"]) | (data_initial["MMFF_InChI"] == data_initial["QM_InChI"])]
            data.to_csv(os.path.join(self.outdir, "data_consistent_loose.csv"), index = False)

    def check_others(self, infile):
        """
        Check radicals and partial charges, MMFF_SMILES is used to check
        
        Args:
            infile: data index file
        """
        update_list = []
        os.chdir(self.datadir)
        data = pd.read_csv(os.path.join(self.outdir, infile))
        index_list = list(data["index"])
        for i in index_list:
            mol = Chem.SDMolSupplier(str(i) + self.suffix)[0]
            if Descriptors.NumRadicalElectrons(mol) == 0:
                if len([atom for atom in mol.GetAtoms() if atom.GetFormalCharge() != 0]) == 0:
                    update_list.append(i)
        data = data[data["index"].isin(update_list)]
        data.to_csv(os.path.join(self.outdir, infile.split(".")[0] + "_rmrpc.csv"), index = False)
        os.chdir(self.olddir)

    
