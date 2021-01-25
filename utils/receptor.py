import mdtraj as mt
import os
import numpy as np
import pandas as pd
from .atomtype import *


def residue_mapper():

    filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../data/non-standard_residues.dat")
    df = pd.read_csv(filename, header=0, sep="\t")
    mapper = {}
    for i in range(df.shape[0]):
        mapper[df.values[i][1]] = df.values[i][1]

    mapper['CYSS'] = 'CYS'

    return mapper


class ProteinParser(object):
    """Featurization of Protein-Ligand Complex based on
    onion-shape distance counts of atom-types.
    Parameters
    ----------
    pdb_fn : str
        The input pdb file name. The file must be in PDB format.
    Attributes
    ----------
    pdb : mdtraj.Trajectory
        The mdtraj.trajectory object containing the pdb.
    receptor_indices : np.ndarray
        The receptor (protein) atom indices in mdtraj.Trajectory
    rec_ele : np.ndarray
        The element types of each of the atoms in the receptor
    pdb_parsed_ : bool
        Whether the pdb file has been parsed.
    distance_computed : bool
        Whether the distances between atoms in receptor and ligand has been computed.
    Examples
    --------
    >>> pdb = ProteinParser("input.pdb")
    >>> pdb.parsePDB('protein and chainid 0')
    >>> pdb.coordinates_
    >>> print(pdb.rec_ele)
    """

    def __init__(self, pdb_fn):
        self.pdb = mt.load_pdb(pdb_fn)

        self.receptor_indices = np.array([])
        self.rec_ele = np.array([])

        self.pdb_parsed_ = False
        self.coordinates_ = None

        self.residue_names = None

    def get_coordinates(self):
        """
        Get the coordinates in the pdb file given the atom indices.
        Returns
        -------
        self : an instance of itself
        """
        # unit: angstrom
        self.coordinates_ = self.pdb.xyz[0][self.receptor_indices] * 10.0 # bug fixed for atomic unit

        return self

    def _residue_standard(self, resname, mapper):
        if resname in all_residues:
            return resname
        else:
            if resname in mapper.keys():
                return mapper[resname]
            else:
                return "OTH"


    def parsePDB(self, rec_sele="protein"):
        """
        Parse the pdb file and get the detail information of the protein.
        Parameters
        ----------
        rec_sele : str,
            The string for protein selection. Please refer to the following link.
        References
        ----------
        Mdtraj atom selection language: http://mdtraj.org/development/atom_selection.html
        Returns
        -------
        """ 

        mapper = residue_mapper()

        # obtain the atom indices of the protein (only protein atoms are selected,
        # water and ions are discarded in this step)
        top = self.pdb.topology
        self.receptor_indices = top.select(rec_sele)
        self.pdb.atom_slice(self.receptor_indices, inplace=True)
        top = self.pdb.topology

        # get a dataframe containing the atom information.
        _table, _bond = top.to_dataframe()

        # fetch the element type of each one of the protein atom
        self.rec_ele = _table['element'][self.receptor_indices].values
        # fetch the coordinates of each one of the protein atom
        self.get_coordinates()

        self.pdb_parsed_ = True

        self.residue_names = [ "{}_{}_{}".format(self._residue_standard(_table['resName'].values[i], mapper), 
                               _table['resSeq'].values[i], _table['chainID'].values[i]) 
                               for i in range(_table.shape[0])]

        return self
