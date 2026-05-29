import numpy as np
from statistics import NormalDist
from scipy.linalg import solve_banded

def european_put_implicit_fdm(S_max, K, T, r, sigma, M = 500, N = 500):
    dtau = T/N 

    S_partition = np.linspace(0, S_max, M+1)
    tau_partition = np.linspace(0, T, N+1)

    V = np.zeros((N+1, M+1))                            #V[n,i] = option value at time to maturity tau_n := tau_partition[n] and stock price S_i := S_partition[i]
    V[0, :] = np.maximum(K - S_partition, 0)            #initial condition at n = 0 (maturity)
    V[:, 0] = K * np.exp(-r * tau_partition)            #first boundary condition: if price of asset is zero, value of option is the discounted strike price
    V[:, -1] = 0                                        #second boundary condition: we assume S_max is sufficiently large so that the option is essentially worthless if the asset reaches S_max

    i = np.arange(1,M)      #vector of interior 'spatial' indices

    a = -0.5 * dtau * ((sigma**2 * i**2) - (r * i))
    b = 1 + dtau * ((sigma**2 * i**2) + r)
    c = -0.5 * dtau * ((sigma**2 * i**2) + (r * i))

    A_bands = np.zeros((3, M-1))
    A_bands[0, 1:] = c[:-1]     #upper diagonal (with leading zero)
    A_bands[1, :] = b           #main diagonal
    A_bands[2, :-1] = a[1:]     #lower diagonal (with trailing zero)

    for n in range(N):
        rhs = V[n, 1:M].copy()
        rhs[0] -= a[0] * V[n+1, 0]
        rhs[-1] -= c[-1] * V[n+1, M]           #this does nothing since right boundary value is 0
           
        V[n+1, 1:M] = solve_banded((1, 1), A_bands, rhs)    #quicker than np.linalg.solve for tridiagonal matrices
    
    return S_partition, tau_partition, V 


def european_put_closed_form(S, K, T, r, sigma):
    if T == 0: 
        return max(K - S, 0)
    if K == 0:
        return 0
    if S == 0:
        return K * np.exp(-r * T)
    
    N = NormalDist().cdf

    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)

    return K * np.exp(-r * T) * N(-d2) - S * N(-d1)


def european_put_closed_form_sampled(S_max, K, T, r, sigma, M = 500, N = 500):
    S_partition = np.linspace(0, S_max, M+1)
    tau_partition = np.linspace(0, T, N+1)
    V = np.zeros((N+1, M+1))      

    for i in range(N+1):
        for j in range(M+1):
            V[i,j] = european_put_closed_form(S_partition[j], K, tau_partition[i], r, sigma)    

    return S_partition, tau_partition, V      