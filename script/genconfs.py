import os, sys
from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints
from rdkit.ML.Cluster import Butina
from rdkit.Chem import PandasTools
import pandas as pd

def gen_conformers(mol, numConfs=1):
    """Generate conformation with local minimization"""
    
    ### generate conf using ETKDG method ###
    ps = AllChem.ETKDG()
    ps.maxAttempts = 1000
    ps.randomSeed = 1
    ps.pruneRmsThresh = 0.1
    ps.numThreads = 0 
    ids = AllChem.EmbedMultipleConfs(mol, numConfs, ps)
    ### Check MMFF parms ###
    if AllChem.MMFFHasAllMoleculeParams(mol):
        ### MMFF optimize ###
        method = "MMFF"
        for cid in ids:
            _ = AllChem.MMFFOptimizeMolecule(mol, confId = cid)
    else:
        ### UFF optimize ###
        method = "UFF"
        for cid in ids:
            _ = AllChem.UFFOptimizeMolecule(mol, confId = cid)
    return list(ids), method

def cluster_conformers(mol, mode="RMSD", threshold=0.2):
    """
    Cluster conf based on heavy atom rmsd 
    Then Butina is used for clustering
    """
    ### get heavy atom idx ###
    heavyatomidx = []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != 1:
            heavyatomidx.append(a.GetIdx())

    ### align on heavy atom for each pair and get dmat ###
    n = mol.GetNumConformers()
    dmat = []
    for i in range(n):
        for j in range(i):
            dmat.append(Chem.rdMolAlign.AlignMol(mol, mol, i, j, atomMap = [(k, k) for k in heavyatomidx]))
    ### clustering ###
    rms_clusters = Butina.ClusterData(dmat, mol.GetNumConformers(), threshold, isDistData=True, reordering=True)
    
    return rms_clusters

def calc_energy(idx, mol, conformerId, method, minimizeIts=0):
    """
    Set minimizeIts to be 0 to turn off min
    since MMFF opt have been done before
    Here, it is used to get MMFF energy
    """
    if method == "MMFF":
        try:
            mp = AllChem.MMFFGetMoleculeProperties(mol)
            ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=conformerId)
            results = {}
            ### default not minimize,since we already did MMFF optimization, conducting minimization here or not is doesn't matter ###
            if minimizeIts > 0:
                ff.Initialize()
                ff.Minimize(maxIts=minimizeIts)
                ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=conformerId)
            results["energy_abs"] = ff.CalcEnergy()
        except:
            ### for some molecules, such as HF, they can't be minimized ###
            results = {}
            results["energy_abs"] = None
    else:
        try:
            ff = AllChem.UFFGetMoleculeForceField(mol)
            results = {}
            if minimizeIts > 0:
                ff.Initialize()
                ff.Minimize(maxIts=minimizeIts)
                ff = AllChem.UFFGetMoleculeForceField(mol)
            results["energy_abs"] = ff.CalcEnergy()
        except:
            ### for some molecules, such as HF, they can't be minimized ###
            results = {}
            results["energy_abs"] = None
        
    return results
    
def runGenerator(index_list, smiles_list, source_data_name, structure_dir=None, numConfs=300, clusterMethod="RMSD", clusterThreshold=0.2):
    """
    Generate conformation as sdf for all smiles in input
    index_list: the list of initial indexes from source data 
    smiles_list: the rdkit smiles generated using initial source data SMILES
    source_data_name: the source data name, such as ccdc, zinc, and pubchem, which is used to record the initial index in that source data
    structure_dir: the directory used to find the 3D structure file, if it is None, the numConfs should not be None
    numConfs: the number of conformations generated before clustering, defaults to 300; if numConfs=None, we conducted local MMFF minimization for structure
    clusterMethod: the distance calculation method used in clustering, defaults to RMSD
    clusterThreshold: the clustering threshold, defaults to 0.2
    notice: we only conducted one conformation generation for each SMILES, but when numConfs=None, we conducted local optimization for each structure
    return:saved *_confs.sdf for conformation generation or *_min.sdf for local minimization with MMFF; Failed_list
    """
    Failed_list = []
    Used_smiles_list = []
    for idx, smiles in enumerate(smiles_list):
        print(idx, smiles)
       
        if numConfs:
            ### check SMILES to make sure only conduct one time conformation generation for each SMILES ###
            if smiles in Used_smiles_list:
                continue
        
        if structure_dir:
            if str(index_list[idx])  + ".sdf" in os.listdir(structure_dir):
                print("Use 3D Structure as Reference:", str(index_list[idx])  + ".sdf")
                mol = Chem.SDMolSupplier(os.path.join(structure_dir, str(index_list[idx])  + ".sdf"), removeHs=False)[0]
                if mol == None:
                    print("Wring 3D Structure!")
                    continue
            else:
                print("Need 3D Structure File Provided")
        else:
            if numConfs == None:
                print("Can't Conducted Structure Local Minimization without Structure!")
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol != None:
                ### add H to initial SMILES ###
                mol = Chem.AddHs(mol)
            else:
                print("Wrong SMILES!")
                continue

        if mol != None: 
            if numConfs:
                ### generated conformations are saved in *_confors.sdf file ###
                w = Chem.SDWriter(str(index_list[idx]) + "_confors.sdf")
                try:
                    conformerIds, method = gen_conformers(mol, numConfs=numConfs)
                except:
                    ### failed cases have been captured by Failed_list ###
                    ### situtation 1: conformation generation process is failed ### 
                    Failed_list.append(idx)
                    continue
                
                if conformerIds == []:
                    ### situation 2: no conformation has been generated ###
                    Failed_list.append(idx)
                    continue
                ### cluster conformations ###
                rmsClusters = cluster_conformers(mol, clusterMethod, clusterThreshold)
 
                conformerPropsDict = {}
                n = 0
                for clusterId in rmsClusters:
                    n = n + 1
                    ### each cluster, we only keep the centroid ###
                    for conformerId in clusterId[:1]:
                        conformerPropsDict[conformerId] = {}
                        ### structure minimization (optional) and energy calculation ###
                        conformerPropsDict[conformerId]["energy_abs"] = calc_energy(idx, mol, conformerId, method)["energy_abs"]
                        ### situation 3: no minimized energy ###
                        if conformerPropsDict[conformerId]["energy_abs"] == None:
                            Failed_list.append(idx)
                            continue
                        conformerPropsDict[conformerId]["SMILES"] = smiles
                        conformerPropsDict[conformerId]["cluster_no"] = n
                        conformerPropsDict[conformerId][source_data_name + "_id"] = index_list[idx]
                        conformerPropsDict[conformerId]["initial_conformation_id"] = conformerId
                        conformerPropsDict[conformerId]["minimize_method"] = method
                        for key in conformerPropsDict[conformerId].keys():
                            mol.SetProp(key, str(conformerPropsDict[conformerId][key]))
                        w.write(mol, confId = conformerId)
                print("The total number of conformers after clustring: " + str(n))
                ### only append smiles in the Used_smiles_list after the successful conformation generation ###
                Used_smiles_list.append(smiles)
            else:
                ### no conformation generation, just minimization of structure ###
                ### local minimized structure is saved in *_min.sdf ###
                w = Chem.SDWriter(str(index_list[idx]) + "_min.sdf")
                conformerPropsDict = {}
                conformerId = 0 
                method = "MMFF"
                conformerPropsDict[conformerId] = {}
                ### here, we need to conduct minimization using calc_energy function, which is controled by minimizeIts=200 ###
                conformerPropsDict[conformerId]["energy_abs"] = calc_energy(idx, mol, conformerId, method, minimizeIts=200)["energy_abs"]
                if conformerPropsDict[conformerId]["energy_abs"] == None:
                    Failed_list.append(idx)
                    continue
                conformerPropsDict[conformerId]["SMILES"] = smiles
                conformerPropsDict[conformerId][source_data_name + "_id"] = index_list[idx]
                conformerPropsDict[conformerId]["initial_conformation_id"] = conformerId
                conformerPropsDict[conformerId]["minimize_method"] = method
                for key in conformerPropsDict[conformerId].keys():
                    mol.SetProp(key, str(conformerPropsDict[conformerId][key]))
                w.write(mol, confId = conformerId)
                print("Finish Local Minimization with MMFF")
        else:
            print("Wrong Structure!")
        
        w.flush()
        w.close()
    
    return Failed_list

def get_index(infile,smiles_name,index_name):
    """ 
    Get index and SMILES list 
    smiles_name: column name of SMILES
    index_name: column name of index
    """
    infile = pd.read_csv(infile)
    smiles = infile[smiles_name].tolist()
    index = infile[index_name].tolist()
    return smiles, index


if __name__ == "__main__":
    size  = "20"
    if size + "_confs" not in os.listdir("."):
        os.mkdir(size + "_confs")
    os.chdir(size + "_confs")
    index_list = []
    smiles_list = []
    for size in ["20"]:
        smiles, indexs = get_index(size)
        index_list.extend(indexs)
        smiles_list.extend(smiles)
    index_list_redo = []
    smiles_list_redo = []
    problem_list = [int(i.rstrip()) for i in open("../redoconfs.csv")]
    for idx, i in enumerate(index_list):
        if i in problem_list:
            index_list_redo.append(i)
            smiles_list_redo.append(smiles_list[idx])
    print(len(index_list_redo))
    structure_dir = "/Users/jianinglu1/Documents/Python_API_2019/my_code/ccdc_structures"
    runGenerator(index_list_redo, smiles_list_redo, structure_dir=structure_dir, numConfs=1000)


