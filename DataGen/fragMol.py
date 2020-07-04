### This is a script to get fragments using Murcko method ###
import os
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold as Murcko
from rdkit.Chem import AllChem
import numpy as np
import ast

def get_fragments(insmiles):
    """Get core and sidechains using Murcko fragmentation method

    Args:
        insmiles (str): SMILES for molecule

    Returns:
        core: scaffold SMILES
        side: side chain SMILES
    """

    mol_prev = Chem.MolFromSmiles(insmiles)
    
    ### get core using Murcko fragmentation ###
    core = Murcko.MurckoScaffoldSmilesFromSmiles(insmiles)
    if core != "":
        mol_core = Chem.MolFromSmiles(core)
    
        ### get sidechains ###
        mol_side = Chem.rdmolops.DeleteSubstructs(mol_prev, mol_core)
        side = Chem.MolToSmiles(mol_side).split(".")
    else:
        side = [insmiles]
        
    return core, side


def get_number(infile, core_column, side_column):
    """Get occurence number after fragmentation, input should be ... + "_addfrags.csv"

    Args:
        infile (str): input file name
        core_column (str): column name for core
        side_column(str): column name for side

    Returns:
        absolute number (dict): the absolute occurence of the fragments
        molecule number (dict): the number of molecules that have certain fragments (max=1 for each molecule)
    """

    inf = pd.read_csv(infile)
    c_id = inf.columns.get_loc(core_column)
    s_id = inf.columns.get_loc(side_column)
    absolute_number = {frag:0 for frag in []}
    molecule_number = {frag:0 for frag in []}

    for idx, i in enumerate(inf.index):
        if str(inf.iloc[i,c_id]) != "nan" and inf.iloc[i,s_id] != "['']":
            fragments = [inf.iloc[i,c_id]] + ast.literal_eval(inf.iloc[i,s_id])
        elif str(inf.iloc[i,c_id]) == "nan":
            fragments = ast.literal_eval(inf.iloc[i,s_id])
        elif inf.iloc[i,s_id] == "['']":
            fragments = [inf.iloc[i,c_id]]
        
        for f in fragments:
            if f not in absolute_number:
                absolute_number[f] = 1
            else:
                absolute_number[f] += 1
        fragments = set(fragments)
        for f in fragments:
            if f not in molecule_number:
                molecule_number[f] = 1
            else:
                molecule_number[f] += 1

        if idx % 10000 == 0:
            print("Finish", idx)

    return absolute_number, molecule_number

def fragMol(dir, infile, column_name, id=None, outfile=None, cal_occurrence=True):
    """Get fragments of molecules

    Args:
        dir: datadir for the input and output file
        infile (str): input file name, contains the SMILES inform
        column_name (str): column_name for SMILES 
        id (str): defaults to None, can be given as the conlumn name of id if infile has id infor
        outfile (str): defaults to None, outfile will be named as "fragments.csv", can be changed using custome name
        cal_occurrence: whether to calculate occurrence number for fragments, if true, will creat a "fragments_occurance.npz" with frags are the name of fragments and occurance are the occurance numbers

    Returns:
        No return, but will generated two files:
        1. File for fragments 
        2. File for original data plus the fragments information for each molecule. Filename is the infile + "_addfrags.csv"
    """
    dir = os.path.realpath(dir)
    df = pd.read_csv(os.path.join(dir,infile))
    smiles_list = list(df[column_name])
    if id:
        id_list = list(df[id])
    else:
        id_list = list(df.index)
    fragments = []
    if outfile:
        out = open(os.path.join(dir, outfile), "w")
    else:
        out = open(os.path.join(dir, "fragments.csv"), "w")

    out.write("mol_id,frags_id,frags_SMILES\n")
    core_list = []
    side_list = []
    if cal_occurrence:
        occurrence = {}
    f_id = 0
    for idx, smiles in enumerate(smiles_list):
        try:
            core, side = get_fragments(smiles)
        except:
            print(id_list[idx], ": Fragmentation failed")
            core_list.append("Problem")
            side_list.append("Problem")
        else:
            core_list.append(core)
            side_list.append(side)
            frags = set([core] + side)
            for f in frags:
                if f not in fragments and f != "":
                    fragments.append(f)
                    out.write(str(id_list[idx]) + "," + str(f_id) + "," + f + "\n")
                    f_id += 1
                    if cal_occurrence:
                        occurrence.update({f:1})
                elif cal_occurrence and f != "":
                    occurrence[f] += 1
        if idx % 10000 == 0:
            print("Finish", idx)
    out.close()

    if cal_occurrence:
        np.savez(os.path.join(dir, "fragments_occurance.npz"), frags=list(occurrence), occurance=list(occurrence.values()))


    df["core_list"] = core_list
    df["side_list"] = side_list
    df.to_csv(os.path.join(dir, infile.split("/")[-1].split(".")[0] + "_addfrags.csv"), index = None)


if __name__ == "__main__":
    dir = "../test/fragmentation/"
    infile = "test_fragments.csv"
    column_name = "SMILES"
    id = "idx"
    fragMol(dir, infile, column_name, id)
    infile_frag = "../test/fragmentation/test_fragments_addfrags.csv"
    absolute_number, molecule_number = get_number(infile_frag, "core_list", "side_list")
    print(absolute_number)
    print(molecule_number)