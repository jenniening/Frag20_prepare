### Generating EFGs is based on in-house package, which is not available now#
import EFGs
from EFGs import mol2frag
from EFGs import cleavage
import rdkit
import torch
import argparse
from rdkit import Chem
import os
import pandas as pd

def get_EFGs(mol, intype, isomeric=True):
    """Get EFGs from mol, mol should be RDKit mol file

    Args:
        mol (molecule object): input molecule
        intype (str): sdf or smiles
        isomeric (boolean): whether to consider isomeric infor. Defaults to True.

    Returns:
        EFGs: Extended functional group for molecules
    """

    if intype == "sdf":
        a,b = mol2frag(mol, TreatHs='include', returnidx=False, toEnd=False, isomericSmiles=isomeric)
    elif intype == "smiles":
        a,b = mol2frag(mol, returnidx=False, toEnd=False, isomericSmiles=isomeric)
    EFGs = set(a + b)
    return EFGs

def get_EFGs_dic(dataset, intype, isomeric=True):
    """Get EFGs_dic, which uses EFGs as keys, and list of index of mols that have this EFG as value
    """

    EFGs_dic = {}
    for idx, mol in enumerate(dataset):
        if intype == "smiles":
            mol = Chem.MolFromSmiles(mol)
        try:
            EFGs = get_EFGs(mol, intype, isomeric=isomeric)
        except:
            print(idx, mol)
            print(Chem.MolToSmiles(mol))
        for i in EFGs:
            if i not in EFGs_dic:
                EFGs_dic[i] = [idx]
            else:
                EFGs_dic[i].append(idx)
        if idx % 10000 == 0:
            print(idx)
    return EFGs_dic

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate EFGs dictionary for dataset')
    parser.add_argument('--dataset', type=str, help="dataset directory")
    #parser.add_argument('--isomeric', type=bool, help="isomeric EFGs")
    parser.add_argument('--outfile', type=str, help="outfile prefix")
    parser.add_argument('--cutoff', type=float, help="cutoff for EFGs")
    args = parser.parse_args()

    dataset = args.dataset
    if dataset.split(".")[-1] == "pt":
        dataset = torch.load(dataset)
        intype = "sdf"

    elif dataset.split(".")[-1] == "csv":
        dataset = pd.read_csv(dataset)
        dataset = list(dataset["smiles_preprocessed"])
        intype = "smiles"

    outfile = args.outfile
    if "add3D" in outfile:
        isomeric=True
    else:
        isomeric=False
    #isomeric = args.isomeric
    print(isomeric)
    cutoff = args.cutoff
    print(cutoff)
    if cutoff == 0:
        EFGs_dic = get_EFGs_dic(dataset, intype=intype, isomeric=isomeric)
        EFGs_dic_frequency = {key:len(value) for key, value in EFGs_dic.items()}
        torch.save(EFGs_dic, outfile + "_EFGs_full.pt")
        torch.save(EFGs_dic_frequency, outfile + "_EFGs_frequency.pt")
    else:
        if outfile + "_EFGs_frequency.pt" in os.listdir("."):
            EFGs_dic_frequency = torch.load(outfile + "_EFGs_frequency.pt")
            print("full EFGs:", len(EFGs_dic_frequency))
            ### get EFGs with certain cutoff ###
            cleavage(EFGs_dic_frequency, alpha=cutoff, isomericSmiles=isomeric)
            print("EFGs after cut:", len(EFGs_dic_frequency))
            torch.save(EFGs_dic_frequency, outfile + "_cutoff" + str(cutoff) + "_EFGs_frequency.pt")
        else:
            EFGs_dic = get_EFGs_dic(dataset, intype=intype, isomeric=isomeric)
            EFGs_dic_frequency = {key:len(value) for key, value in EFGs_dic.items()}
            torch.save(EFGs_dic, outfile + "_EFGs_full.pt")
            torch.save(EFGs_dic_frequency, outfile + "_EFGs_frequency.pt")
            EFGs_dic_frequency = torch.load(outfile + "_EFGs_frequency.pt")
            print("full EFGs:", len(EFGs_dic_frequency))
            cleavage(EFGs_dic_frequency, alpha=cutoff, isomericSmiles=isomeric)
            print("EFGs after cut:", len(EFGs_dic_frequency))
            torch.save(EFGs_dic_frequency, outfile + "_cutoff" + str(cutoff) + "_EFGs_frequency.pt")










