# Frag20_prepare
This is a package for data generation including useful functions for **molecule fragmentation**, **molecule selection**, **1D to 3D labelling**: conformation generation, file conversion, MMFF and QM calculation, and **dataset preparation**

### Prerequisites
Python <br>
Numpy (used to save data) <br>
Pandas (used to read data) <br>
Openbabel (used to convert file into different format)<br>
```
conda install -c openbabel openbabel
```

RDkit (used to generate conformations)<br>
```
conda install -c rdkit rdkit
```

Pytorch (used to save data) <br>
Install follow the [link](https://pytorch.org/?utm_source=Google&utm_medium=PaidSearch&utm_campaign=%2A%2ALP+-+TM+-+General+-+HV+-+US&utm_adgroup=Install+PyTorch&utm_keyword=install%20pytorch&utm_offering=AI&utm_Product=PyTorch&gclid=EAIaIQobChMIoJbb_Nmy6gIViorICh0PigPMEAAYASAAEgLgy_D_BwE) <br><br>
Ase (used to convert unit) <br> 
```
conda install -c conda-forge ase
```

EFGs (used to generated extend functional groups, developed by Jieyu Lu)<br>

Please check the link for detailed information: https://github.com/HelloJocelynLu/EFGs

### Setup
Install our DataGen package
```
python setup.py install
```

### Tutorials
#### 1. [Molecule fragmentation and selection](http://htmlpreview.github.io/?https://github.com/jenniening/Frag20_prepare/blob/master/tutorials/Data_preparation_tutorial.html)
This is tutorial for basic functions we used to cut molecules into fragments, generate EFGs, and select molecules.
#### 2. [1D to 3D labelling and dataset prepartion](http://htmlpreview.github.io/?https://github.com/jenniening/Frag20_prepare/blob/master/tutorials/1D_to_3D_tutorial.html)
This tutorial includes all steps in 1D to 3D labelling and also the functions used for prepare input dataset for PhysNet model.

### Reference
Jianing Lu, Song Xia, Jieyu Lu and Yingkai Zhang. *J. Chem. Inf. Model.*, **61**. (2021)
https://pubs.acs.org/doi/10.1021/acs.jcim.1c00007


