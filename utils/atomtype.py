
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