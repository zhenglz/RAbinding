
import numpy as np
import mdtraj as mt
import itertools


def fast_distance_pairs(coords_pro, coords_lig):
    return np.array([np.linalg.norm(coord_pro - coords_lig, axis=-1) for coord_pro in coords_pro]).ravel()


def distance_pairs_mdtraj(coord_pro, coord_lig):
    xyz = np.concatenate((coord_pro, coord_lig), axis=0)
    #print(xyz.shape)
    # mdtraj use nanometer for coordinations,
    # convert angstrom to nanometer
    xyz = xyz.reshape((1, -1, 3)) * 0.1
    # for the xyz, the sample is (N_frames, N_atoms, N_dim)
    # N_frames, it is usually 1 for a normal single-molecule PDB
    # N_atoms is the number of atoms in the pdb file
    # N_dims is the dimension of coordinates, 3 here for x, y and z

    # create a mdtraj Trajectory object,Topology object could be ignored.
    traj = mt.Trajectory(xyz=xyz, topology=None)
    # create a list of atom-pairs from atom index of protein and ligand
    atom_pairs = itertools.product(np.arange(coord_pro.shape[0]),
                                   np.arange(coord_pro.shape[0], coord_pro.shape[0]+coord_lig.shape[0]))

    # convert the distance to angstrom from nanometer.
    # Actually we could just leave it as angstrom from the beginning for faster calculation,
    # but it is more reasonable to do it in order to aligning with mdtraj-style calculation.
    dist = mt.compute_distances(traj, atom_pairs)[0] * 10.0

    return dist


def residue_min_distance(coord_pro, coord_lig, residue_names, receptor_elements, ligand_elements):
    # combine same residues in the residue_names list
    uniq_res = []
    for _residue in residue_names:
        if _residue not in uniq_res:
            uniq_res.append(_residue)
    
    print(coord_pro.shape, len(residue_names))

    #assert coord_pro.shape[0] == len(residue_names)
    assert coord_pro.shape[0] == len(receptor_elements)
    assert coord_lig.shape[0] == len(ligand_elements)

    _ligand_indices = np.array([x for x in range(len(ligand_elements)) if ligand_elements[x] != "H"])

    results = np.zeros((len(uniq_res), _ligand_indices.shape[0]))

    for i, resid in enumerate(uniq_res):
        _receptor_indices = np.array([x for x in range(coord_pro.shape[0]) 
                                      if (residue_names[x] == resid and 
                                          receptor_elements[x] != "H")])
        #print(_receptor_indices, _ligand_indices)
        
        _distances = fast_distance_pairs(coord_pro[_receptor_indices], coord_lig[_ligand_indices])
        _distances = _distances.reshape((-1, _ligand_indices.shape[0]))
        
        _min_dist = np.mean(_distances, axis=0)
        #print(_min_dist)
        results[i] = _min_dist
    
    return (results, uniq_res, [x for x in ligand_elements if x != "H"])



