import numpy as np 

def compute_greeks(S_partition, V):
    dS = S_partition[1] - S_partition[0]

    Delta = np.full_like(V, np.nan, dtype = float)
    Gamma = np.full_like(V, np.nan, dtype = float)

    Delta[:, 1:-1] = (V[:, 2:] - V[:, :-2]) / (2 * dS)
    Gamma[:, 1:-1] = (V[:, 2:] - 2*V[:, 1:-1] + V[:, :-2]) / (dS ** 2)

    return Delta, Gamma 