#!/usr/bin/ python 3.X

import time
import numpy as np
import pandas as pd
from argparse import RawDescriptionHelpFormatter
import argparse
import mdtraj as md
import itertools
from collections import OrderedDict
from utils.atomtype import *
from utils.receptor import ProteinParser
from utils.ligand import LigandParser
from utils.distances import residue_min_distance
from utils.contacts import dist2count_simple
from utils.parallel import ParallelSim


# all residue-atom combination pairs
keys = ["_".join(x) for x in list(itertools.product(all_residues, all_elements))]


def generate_features(receptor_file, ligand_file, cutoffs):
    # process ligand
    ligand = LigandParser(ligand_file)
    ligand.parseMol2()

    # process receptor
    receptor = ProteinParser(receptor_file)
    receptor.parsePDB()

    # get receptor information
    _residue_names = receptor.residue_names
    #print(_residue_names)
    _residue_elements = receptor.rec_ele
    _residue_coordinates = receptor.coordinates_
    #print("RECEPTOR COORD", _residue_coordinates.shape)

    # get ligand information
    _ligand_elements = ligand.lig_ele
    #print(_ligand_elements)
    _ligand_coordinates = ligand.coordinates_
    #print("LIGAND COORD", _ligand_coordinates.shape)

    # get mini distances
    _min_distances, _residues, _non_hydrogen_elements = residue_min_distance(_residue_coordinates, _ligand_coordinates, 
                                                                      _residue_names, _residue_elements, _ligand_elements)

    # Types of the residue and the atom 
    #new_residue = list(map(get_residue, cplx.residues))
    new_residue = [x.split("_")[0] for x in _residues]
    new_lig = list(map(get_elementtype, _non_hydrogen_elements))

    # residue-atom pairs
    residues_lig_atoms_combines = ["_".join(x) for x in list(itertools.product(all_residues, all_elements))]
    #print("LENGTH RESIDUE COMBINES ", len(residues_lig_atoms_combines), len(set(residues_lig_atoms_combines)))
    #print(residues_lig_atoms_combines)

    # calculate the number of contacts in different shells
    counts = []
    onion_counts = []
    final_results = OrderedDict()

    for i, cutoff in enumerate(cutoffs):
        counts_ = dist2count_simple(_min_distances, cutoff)#.ravel()
        counts.append(counts_.ravel())
        #print(len(onion_counts), counts_.shape, )
        if i == 0:
            onion_counts.append(list(counts_.ravel()))
        else:
            onion_counts.append(list(counts_.ravel() - counts[-1]))

        for j, _key in enumerate(residues_lig_atoms_combines):
            #print("COMBINE SIZE", len(residues_lig_atoms_combines))
            _new_key = _key + "_" + str(i)
            final_results[_new_key] = 0.0
        
        for j in range(len(new_residue)):
            for k in range(len(new_lig)):
                _new_key = new_residue[j] + "_" + new_lig[k] + "_" + str(i)
                final_results[_new_key] += counts_[j, k]
    
    #print("Number of keys in results ", len(final_results.keys()))

    return list(final_results.values()), receptor_file, ligand_file


if __name__ == "__main__":

    print("Start Now ... ")
    start_time = time.time()

    d = """
        Generate the residue-atom contact features.

        The distance calculated in this script is the minimum distance between residues and atoms,\n
        that is, the distance between the atom (in the ligand) and the closest heavy atom in the specific residue.

        usage:
            python generate_features.py -inp input_complex.dat -out output_features.csv
        """

    parser = argparse.ArgumentParser(description=d, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-inp", type=str, default="input_complex.dat",
                        help="Input. This file includes the path of the pdb file \n"
                             "for the protein-ligand complexes. \n"
                             "Each line in the dat file contains only one pdb file path.")
    parser.add_argument("-out", type=str, default="output_features.csv",
                        help="Output. The output file format is .csv \n"
                             "Each line in the output file contains only the features \n"
                             "of one protein-ligand complex.")
    parser.add_argument("-shells", type=int, default=62,
                        help="Input. The total number of shells. The optional range here \n"
                             "is 10 <= N <= 90.")
    parser.add_argument("-cpu", type=int, default=12,
                        help="Input (Optional). Number of CPUs for parallel computation.")                        

    args = parser.parse_args()

    with open(args.inp) as lines:
        inputs = [x.split()[:2] for x in lines if ("#" not in x and len(x))]

    index = []
    outermost = 0.5 * (args.shells + 1)
    ncutoffs = np.linspace(2.0, outermost, args.shells)

    l = len(inputs)
    if args.cpu <= 1:
        results = []
        for i, fns in enumerate(inputs):
            rec, lig = fns
            print(rec, lig)
            result, _r, _l = generate_features(rec, lig, ncutoffs)
            results.append(list(np.array(result).ravel()))
            index.append(_r+ "_"+ _l)
            print("INFO: data shape = {} | index {} | total {}".format(np.array(results).shape, i, l))
    else:
        paral = ParallelSim(args.cpu)
        #arguments = []
        results = []
        for i, fns in enumerate(inputs):
            rec, lig = fns
            _args = [rec, lig, ncutoffs]
            paral.add(generate_features, _args)
        paral.run()
        all_results = paral.get_results()
        results = [x[0] for x in all_results]
        index = [x[1] + "_" + x[2] for x in all_results]
        print("INFO: data shape = {} | index {} | total {}".format(np.array(results).shape, i, l))

    columns = []
    for i, n in enumerate(keys * len(ncutoffs)):
        columns.append(n + '_' + str(i))

    df = pd.DataFrame(np.array(results), index=index, columns=columns)
    df.to_csv(args.out, float_format='%.1f')

    end_time = time.time()
    print("Total computation time: %.3f seconds" % (end_time - start_time))
