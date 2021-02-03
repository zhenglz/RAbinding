
elements_ligand = ["H", "C", "O", "N", "P", "S", "Hal", "DU"]
elements_protein = ["H", "C", "O", "N", "P", "S", "Hal", "DU"]


# Define all residue types
all_residues = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'TYR', 'TRP', 'SER',
               'THR', 'CYS', 'MET', 'ASN', 'GLN', 'ASP', 'GLU', 'LYS', 'ARG', 'HIS', 'OTH']
               
def get_residue(residue):
    if residue in all_residues:
        return residue
    else:
        return 'OTH'

# Define all element types
all_elements = ['H', 'C',  'O', 'N', 'P', 'S', 'Hal', 'DU']
Hal = ['F', 'Cl', 'Br', 'I']

def get_elementtype(e):

    if e in all_elements:
        return e
    elif e in Hal:
        return 'Hal'
    else:
        return 'DU'


def get_protein_elementtype(e):
    if e in elements_protein:
        return e
    else:
        return "DU"


def get_ligand_elementtype(e):
    '''if e == "C.ar":
        return "CAR"  '''
    if e.split(".")[0] in elements_ligand:
        return e.split(".")[0]
    else:
        return "DU"


class Molecule(object):
    """Small molecule parser object with Rdkit package.
    Parameters
    ----------
    in_format : str, default = 'smile'
        Input information (file) format.
        Options: smile, pdb, sdf, mol2, mol
    Attributes
    ----------
    molecule_ : rdkit.Chem.Molecule object
    mol_file : str
        The input file name or Smile string
    converter_ : dict, dict of rdkit.Chem.MolFrom** methods
        The file loading method dictionary. The keys are:
        pdb, sdf, mol2, mol, smile
    """

    def __init__(self, in_format="smile"):

        self.format = in_format
        self.molecule_ = None
        self.mol_file = None
        self.converter_ = None
        self.mol_converter()

    def mol_converter(self):
        """The converter methods are stored in a dictionary.
        Returns
        -------
        self : return an instance of itself
        """
        self.converter_ = {
            "pdb": Chem.MolFromPDBFile,
            "mol2": Chem.MolFromMol2File,
            "mol": Chem.MolFromMolFile,
            "smile": Chem.MolFromSmiles,
            "sdf": Chem.MolFromMolBlock,
            "pdbqt": self.babel_converter,
        }

        return self

    def babel_converter(self, mol_file, output):
        if os.path.exists(mol_file):
            try:
                cmd = 'obabel %s -O %s > /dev/null' % (mol_file, output)
                job = sp.Popen(cmd, shell=True)
                job.communicate()

                self.molecule_ = self.converter_['pdb']()
                return self.molecule_
            except:
                return None
        else:
            print("No such input for converting: ", mol_file)

        return None

    def load_molecule(self, mol_file):
        """Load a molecule to have a rdkit.Chem.Molecule object
        Parameters
        ----------
        mol_file : str
            The input file name or SMILE string
        Returns
        -------
        molecule : rdkit.Chem.Molecule object
            The molecule object
        """

        self.mol_file = mol_file
        if not os.path.exists(self.mol_file):
            print("Molecule file not exists. ")
            return None

        if self.format not in ["mol2", "mol", "pdb", "sdf", "pdbqt"]:
            print("File format is not correct. ")
            return None
        else:
            try:
                self.molecule_ = self.converter_[self.format](self.mol_file, sanitize=True,)
            except RuntimeError:
                return None

            return self.molecule_
