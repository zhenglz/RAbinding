import mdtraj as mt
from biopandas.mol2 import PandasMol2
import os
import numpy as np
import pandas as pd
from .atomtype import *


class LigandParser(object):
    """Parse the ligand with biopanda to obtain coordinates and elements.
    Parameters
    ----------
    ligand_fn : str,
        The input ligand file name.
    Methods
    -------
    Attributes
    ----------
    lig : a biopandas mol2 read object
    lig_data : a panda dataframe object holding the atom information
    coordinates : np.ndarray, shape = [ N, 3]
        The coordinates of the atoms in the ligand, N is the number of atoms.
    """

    def __init__(self, ligand_fn):
        self.lig_file = ligand_fn
        self.lig = None
        self.lig_data = None

        self.lig_ele = None
        self.coordinates_ = None
        self.ligand_parsed_ = False

        self.ligand_format = os.path.basename(ligand_fn).split(".")[-1]

    def _format_convert(self, input, output):
        mol = Molecule(in_format=input.split(".")[-1])
        mol.babel_converter(input, output)
        return self

    def get_element(self):
        ele = list(self.lig_data["atom_type"].values)
        self.lig_ele = list(map(get_ligand_elementtype, ele))
        return self

    def get_coordinates(self):
        """
        Get the coordinates in the pdb file given the ligand indices.
        Returns
        -------
        self : an instance of itself
        """
        self.coordinates_ = self.lig_data[['x', 'y', 'z']].values
        return self

    def parseMol2(self):
        try:
            self.lig = PandasMol2().read_mol2(self.lig_file)
        except ValueError:
            print("INFO: Warning, parse mol2 file error, converting to PDB instead ......")
            templ_ligfile = self.lig_file + "templ.pdb"
            # convert mol2 format to pdb format with rdkit
            self._format_convert(self.lig_file, templ_ligfile)

            if os.path.exists(templ_ligfile):
                self.parsePDB(templ_ligfile)
                os.remove(templ_ligfile)
                return self

        self.lig_data = self.lig.df
        self.get_element()
        self.get_coordinates()
        self.ligand_parsed_ = True

        return self
    
    def parsePDB(self):
        self.lig = mt.load_pdb(self.lig_file)
        top = self.lig.topology
        table, bond = top.to_dataframe()
        # get the element-type of all atoms
        self.lig_ele = [x if x in elements_ligand else "DU" for x in list(table['element'])]

        # mdtraj use nanometer as coordinates unit.
        # nano-meter to angstrom
        self.coordinates_ = self.lig.xyz[0] * 10.0
        self.lig_data = table
        self.lig_data['x'] = self.coordinates_[:, 0]
        self.lig_data['y'] = self.coordinates_[:, 1]
        self.lig_data['z'] = self.coordinates_[:, 2]
        self.ligand_parsed_ = True

        return self

    def parsePDBQT(self):
        self.lig = PdbqtParser(self.lig_file)
        #self.lig._get_atom_lines()

        self.lig_data = self.lig.to_dataframe()
        self.lig_ele = self.lig_data['element']
        self.coordinates_ = np.zeros((self.lig_data.shape[0], 3))
        self.coordinates_[:, 0] = self.lig_data['x']
        self.coordinates_[:, 1] = self.lig_data['y']
        self.coordinates_[:, 2] = self.lig_data['z']

        return self


class PdbqtParser(object):

    """Parse PDBQT file.
    Parameters
    ----------
    pdbqt_fn : str,
        Input pdbqt file.
    Examples
    --------
    >>> pdbqt = PdbqtParser("pdb.pdbqt")
    >>> df = pdbqt.to_dataframe()
    >>> df.values
    >>> df.columns
    >>> df['serial']
    """

    def __init__(self, pdbqt_fn=""):
        self.fn = pdbqt_fn

    def _get_atom_lines(self):
        if not os.path.exists(self.fn):
            print("INFO: No such pdbqt ")
            return []

        with open(self.fn) as lines:
            lines = [x for x in lines if x.startswith("ATOM")]

        return lines

    def _element_determinator(self, element, atomname):

        if len(element) == 1:
            if element == "A":
                return "CAR"
            elif element.upper() not in ['C', 'H', 'N', 'O', 'P', 'S', 'K', 'I', 'F']:
                print("INFO: find unusal element ", element, "and it will be replace by %s" % atomname[0].upper())
                return atomname[0].upper()
            else:
                return element.upper()
        elif len(element) == 2:
            e = element.upper()
            if e in ['BR', 'CL', 'MG', 'ZN']:
                return e
            else:
                return e[0]
        else:
            return "NA"

    def _parse_atom_line(self, line):
        _atom_id = int(line[6:11].strip())
        _atom_name = line[12:16].strip()
        _chainid = line[21]
        _res_name = line[17:20].strip()
        try:    
            _res_id = int(line[22:26].strip())
        except ValueError:
            _res_id = 0
        # coordination in unit: angstrom
        _x = float(line[30:38].strip())
        _y = float(line[38:46].strip())
        _z = float(line[46:54].strip())
        _element = self._element_determinator(line[76:79].strip(), _atom_name)
        try:
            _charge = float(line[70:76].strip())
        except ValueError:
            _charge = 0.0        

        return [_atom_id, _atom_name, _chainid, _res_name,
                _res_id, _x, _y, _z, _element, _charge]

    def to_dataframe(self):
        atom_lines = self._get_atom_lines()
        if atom_lines is None:
            return None
        else:
            _data = []
            for line in atom_lines:
                _data.append(self._parse_atom_line(line))

            _df = pd.DataFrame(_data, columns=['serial', 'name', 'chainID', 'resName',
                                               'resSeq', 'x', 'y', 'z', 'element', 'charge'])

            return _df