
import numpy as np


def distance2counts(megadata):
    d, c, charge_pairs, distance_mode, distance_delta = megadata

    if charge_pairs is None:
        return np.sum((np.array(d) <= c) * 1.0)
    else:
        atompair_in_dist_range = ((np.array(d) <= c) & (np.array(d) < c - distance_delta)) * 1.0

        if distance_mode == 'cutoff':
            return np.multiply(atompair_in_dist_range, np.array(charge_pairs) / c)
        else:
            return np.multiply(atompair_in_dist_range, np.divide(charge_pairs, d))

def dist2count_simple(distances, cutoff):
    return (distances <= cutoff) * 1.0

