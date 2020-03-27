import numpy as np
from ase.units import Hartree, eV, Bohr, Ang
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem 
import torch
import os


class data_infor:
    def __init__(self, index, datadir, reference, ftype="confs"):
        """
        Use to get information we need for each molecule
        
        :param index: the index of current structure
        :param datadir: the directory for all structure file and property file
        :param reference: atomic reference energy for QM method
        :param ftype: if local minimization, should be "cry", else, should be "confs", defaults to "confs"

        """
        self._datadir = datadir
        self._index = index
        self._reference = np.load(reference)["atom_ref"]
        ### log file is the Gaussin output file ###
        self._log = os.path.join(self._datadir, index + ".opt.log")
        self._QMsdf = os.path.join(self._datadir, index + ".opt.sdf")
        if ftype == "cry":
            self._MMFFsdf = os.path.join(self._datadir, index + "_min.sdf")
        else:
            self._MMFFsdf = os.path.join(self._datadir, index + ".sdf")
        self._loglines = open(self._log).readlines()
        self._QMlines = open(self._QMsdf).readlines()
        self._MMFFlines = open(self._MMFFsdf).readlines()
        self._QMmol = Chem.SDMolSupplier(self._QMsdf)[0]
        self._MMFFmol = Chem.SDMolSupplier(self._MMFFsdf)[0]
        ### elementdict is for conversion of element, which is the atomic number of elemnet ###
        self._elementdict = {'H':1, 'B':5, 'C':6, 'N':7, 'O':8, 'F':9, 'P':15, 'S':16, 'Cl':17, 'Br':35}
        ### elementdict_periodic is for conversion of element_p, which is the column and period of element ###
        self._elementdict_periodic = {"H":[1,1], "B":[2,3], "C":[2,4], "N":[2,5], "O":[2,6], "F":[2,7], "P":[3,5], "S":[3,6], "Cl":[3,7], "Br":[4,7]}
        ### initial target names we extracted from Gaussian output file ###
        self._target_names = ['A', 'B', 'C', 'mu', 'alpha', 'ehomo', 'elumo', 'egap', 'R2', 'zpve', 'U0', 'U', 'H', 'G', 'Cv', 'E']
        ### conversions is for unit conversion, here, we converted all energy unit from Hartree into eV, which is same as the unit in atom reference energy###
        self._conversions = [ 1., 1., 1., 1., 1.,Hartree/eV, Hartree/eV, Hartree/eV, 1., Hartree/eV, Hartree/eV, Hartree/eV, Hartree/eV, Hartree/eV, 1., Hartree/eV]
        ### init_target is the target directly extracted from Gaussian log file ###
        self._init_target = None
        ### target is the target after conversion and last four are atomiztion energies for U0, U, H, G ###
        self._target = None
        ### number of atoms ###
        self._natoms = None
        ### element atomic number ###
        self._elements = None
        ### element with column and period ###
        self._elements_p = None
        ### mulliken charge ###
        self._charges_mulliken = None
        ### coords from MMFF optimized geomerty ###
        self._MMFFcoords = None
        ### coords from QM optimized geometry ###
        self._QMcoords = None
        ### molecule dipole vector [x, y, z] ###
        self._dipole = None

    @property
    def datadir(self):
        return self._datadir

    @property
    def index(self):
        return self._index

    @property
    def reference(self):
        return self._reference

    @property
    def elementdict(self):
        return self._elementdict

    @property
    def elementdict_periodic(self):
        return self._elementdict_periodic

    @property
    def target_names(self):
        return self._target_names

    @property
    def conversions(self):
        return self._conversions

    @property
    def log(self):
        return self._log

    @property 
    def QMsdf(self):
        return self._QMsdf

    @property 
    def MMFFsdf(self):
        return self._MMFFsdf

    @property 
    def loglines(self):
        return self._loglines
    
    @property
    def QMlines(self):
        return self._QMlines

    @property
    def MMFFlines(self):
        return self._MMFFlines
    @property 
    def QMmol(self):
        return self._QMmol

    @property 
    def MMFFmol(self):
        return self._MMFFmol

    @property 
    def init_target(self):
        if self._init_target == None:
            self.get_initial_target()
        return self._init_target
    @property
    def target(self):
        if self._target == None:
            self.get_target()
        return self._target

    @property
    def natoms(self):
        if self._natoms == None:
            self.get_elements()
        return self._natoms

    @property
    def elements(self):
        if self._elements == None:
            self.get_elements()
        return self._elements

    @property
    def elements_p(self):
        if self._elements_p == None:
            self.get_elements()
        return self._elements_p

    @property
    def charges_mulliken(self):
        if self._charges_mulliken == None:
            self.get_mulliken_charges()
        return self._charges_mulliken

    @property
    def MMFFcoords(self):
        if self._MMFFcoords == None:
            self.get_coordinates()
        return self._MMFFcoords

    @property
    def QMcoords(self):
        if self._QMcoords == None:
            self.get_coordinates()
        return self._QMcoords

    @property
    def dipole(self):
        if self._dipole == None:
            self.get_dipole()
        return self._dipole

    def get_initial_target(self):
        """ 
        Get the information from Gaussian 09 output xx.log, no conversion
        Properties list contains 16 values. 
        The first 15 is the same as in QM9 and the last value is E(B3LYP) (Ha unit)
        """
        log = self._loglines
        propdict = {i:None for i in self._target_names}
        for idx, item in enumerate(log):
            if item.startswith(' Rotational constants'):
                vals = item.split()
                propdict['A'] = vals[-3]
                propdict['B'] = vals[-2]
                propdict['C'] = vals[-1]
            elif item.startswith(' Dipole moment'):
                propdict['mu'] = log[idx+1].split()[-1]
            elif item.startswith(' Isotropic polarizability'):
                propdict['alpha'] = item.split()[-2]
            elif (item.startswith(' Alpha  occ. eigenvalues') and
                log[idx+1].startswith(' Alpha virt. eigenvalues')):
                propdict['ehomo'] = log[idx+1].split()[4]
                propdict['elumo'] = item.split()[-1]
                propdict['egap'] = str(round(float(propdict['ehomo']) - float(propdict['elumo']),5))
            elif item.startswith(' Electronic spatial extent'):
                propdict['R2'] = item.split()[-1]
            elif item.startswith(' Zero-point correction'):
                propdict['zpve'] = item.split()[-2]
            elif item.startswith(' Sum of electronic and zero-point Energies'):
                propdict['U0'] = item.split()[-1]
            elif item.startswith(' Sum of electronic and thermal Energies'):
                propdict['U'] = item.split()[-1]
            elif item.startswith(' Sum of electronic and thermal Enthalpies'):
                propdict['H'] = item.split()[-1]
            elif item.startswith(' Sum of electronic and thermal Free Energies'):
                propdict['G'] = item.split()[-1]
            elif item.startswith(' Total       '):
                propdict['Cv'] = item.split()[-2]
            elif item.startswith(' SCF Done'):
                propdict['E'] = item.split()[4]
        try:
            propvals = [float(propdict[i]) for i in self._target_names]
        except:
            propvals = []
            for i in self._target_names:
                if propdict[i] != None and propdict[i] != "************" and propdict[i] != "2.846627976811D+03":
                    propvals.append(float(propdict[i]))
                else:
                    propvals.append(None)
        self._init_target = propvals

    def get_target(self):
        """ Get properties for each molecule, and convert properties in Hartree unit into eV unit """
        conversions_dict = {}
        for i in range(16):
            conversions_dict[i] = self._conversions[i]
        self._target = [float(self.init_target[i])*conversions_dict[i] for i in range(16)]
        reference_total_1 = np.sum([self._reference[i][1] for i in self.elements])
        reference_total_2 = np.sum([self._reference[i][2] for i in self.elements])
        reference_total_3 = np.sum([self._reference[i][3] for i in self.elements])
        reference_total_4 = np.sum([self._reference[i][4] for i in self.elements])
        self._target.append(self._target[10] - reference_total_1)
        self._target.append(self._target[11] - reference_total_2)
        self._target.append(self._target[12] - reference_total_3)
        self._target.append(self._target[13] - reference_total_4)

    def get_elements(self):
        """ Get elements infor for both atomic number, and periodic based """
        QMnatoms = int(self._QMlines[3].split()[0])
        MMFFnatoms = int(self._MMFFlines[3].split()[0])
        assert QMnatoms == MMFFnatoms, "Error: different number of atoms in mmff and qm optimized files"
        self._natoms = QMnatoms
        atoms = self._QMlines[4:self._natoms + 4]
        elements = []
        elements_p = []
        for atom in atoms:
            atom = atom.split()
            elements.append(self._elementdict[atom[3]])
            elements_p.append(self._elementdict_periodic[atom[3]])
        self._elements = elements
        self._elements_p = elements_p

    def get_coordinates(self):
        """ Get atom coordinates for both MMFF and QM """
        atoms_MMFF = self._MMFFlines[4:self.natoms + 4]
        atoms_QM = self._QMlines[4:self.natoms + 4]
        positions_QM = []
        positions_MMFF = []
        for atom in atoms_QM:
            atom = atom.split()
            positions_QM.append([float(pos) for pos in atom[:3]])

        for atom in atoms_MMFF:
            atom = atom.split()
            positions_MMFF.append([float(pos) for pos in atom[:3]])
        
        self._MMFFcoords = positions_MMFF
        self._QMcoords = positions_QM
        
    def get_mulliken_charges(self):
        """ Get Mulliken charges """
        index = [ idx for idx, line in enumerate(self._loglines) if line.startswith(" Mulliken charges:")][0] + 2
        natoms_old = self.natoms
        natoms = self.natoms
        try: 
            charges = [float(line.split()[-1]) for line in self._loglines[index: index + natoms]]
        except:
            charges = []
            for idx, line in enumerate(self._loglines[index:]):
                if idx < natoms:
                    ### remove calculation comments in Mulliken charges part ###
                    try: 
                        charge = float(line.split()[-1])
                        charges.append(charge) 
                    except:
                        print(line)
                        natoms += 1
                        continue
        assert len(charges) == natoms_old, "Error: charges are wrong"
        self._charges_mulliken = charges
    
    def get_dipole(self):
        """ Calculate dipole using coordinates and charge for each atom """
        coords = self.QMcoords
        dipole = [[coords[i][0] * self.charges_mulliken[i],coords[i][1] * self.charges_mulliken[i], coords[i][2] * self.charges_mulliken[i]] for i in range(self.natoms)]
        dipole = np.sum(dipole,axis = 0)
        self._dipole = dipole  
    
            
def prepare_PhysNet_input(index_list, output, datadir, reference, method="QM", ftype="confs", largest_num_atoms=29):
    """
    Generate standard PhysNet input numpy file 
    
    :param index_list: the list for all input indexes
    :param output: output name for the numpy file
    :param datadir: directory for all sdf and log files
    :param reference: atomic reference energies for QM method, which is used to get the atomization energies
    :param method: if use QM geometries, method = "QM"; use MMFF geometries, method = "MMFF"
    :param ftype: if local minimization, should be "cry", else, should be "confs", defaults to "confs"
    :param largest_num_atoms: the largest number of atoms of dataset, defaults to 29(QM9), should be got before the data prepartion.
    notice: this is for neutral, equilibrium molecules (at local minimization point).
    """
    R_list,Q_list,D_list,E_list,F_list,Z_list,N_list = [], [], [], [], [], [], []
    for idx, i in enumerate(index_list):
        coords = [[0.0,0.0,0.0] for i in range(largest_num_atoms)]
        atoms = [0 for i in range(largest_num_atoms)]
        total_charge = 0
        tmp_data = data_infor(i, datadir, reference, ftype)
        num_atoms = tmp_data.natoms
        if method == "QM":
            coords_new = tmp_data.QMcoords
        else:
            coords_new = tmp_data.MMFFcoords
        coords[0:num_atoms] = coords_new
        dipole = tmp_data.dipole
        atoms_new = tmp_data.elements
        atoms[0:num_atoms] = atoms_new
        energy = tmp_data.target[16]
        force = [[0,0,0] for i in range(largest_num_atoms)]
    
        R_list.append(coords)
        Q_list.append(total_charge)
        D_list.append(dipole)
        E_list.append(energy)
        F_list.append(force)
        Z_list.append(atoms)
        N_list.append(num_atoms)
        if idx % 10000 == 0 and idx != 0:
            print(idx)

    np.savez(output, R=R_list, Q=Q_list, D=D_list, E=E_list, F=F_list, Z=Z_list, N=N_list)

def prepare_torch(index_list, output, datadir, reference, ftype="confs", method="QM"):
    """
    Save rdkit mols into torch file
    
    :param index_list: the list for all input indexes
    :param output: output name for the numpy file
    :param datadir: directory for all sdf and log files
    :param reference: atomic reference energies for QM method, which is used to get the atomization energies
    :param ftype: if local minimization, should be "cry", else, should be "confs", defaults to "confs"
    :param method: optimization method, defaults to "QM"

    """
    sdflist = []
    for idx, i in enumerate(index_list):
        tmp_data = data_infor(i, datadir, reference, ftype)
        if method == "QM":
            sdflist.append(tmp_data.QMmol)
        else:
            sdflist.append(tmp_data.MMFFmol)
        if idx % 10000 == 0 and idx != 0:
            print("Finish ", idx)
    torch.save(sdflist, output)

def prepare_target_csv(index_list, output, datadir, reference, ftype="confs",):
    """
    Save targets into csv file
    
    :param index_list: the list for all input indexes
    :param output: output name for the numpy file
    :param datadir: directory for all sdf and log files
    :param reference: atomic reference energies for QM method, which is used to get the atomization energies
    :param ftype: if local minimization, should be "cry", else, should be "confs", defaults to "confs"

    """
    out = open(output, "w")
    header = ['index','A', 'B', 'C', 'mu', 'alpha', 'ehomo', 'elumo', 'egap', 'R2', 'zpve', 'U0', 'U', 'H', 'G', 'Cv', 'E', 'U0_atom', "U_atom", "H_atom", "G_atom"]
    out.write(",".join(header) + "\n")
    for idx, i in enumerate(index_list):
        tmp_data = data_infor(i, datadir, reference, ftype)
        out.write(str(i) + "," + ",".join([str(i) for i in tmp_data.target]) + "\n")
        if idx % 10000 == 0 and idx != 0:
            print("Finish ", idx) 
    out.close()


