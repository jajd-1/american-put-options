import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
from statistics import NormalDist

def european_put_implicit_fdm(S_max, K, T, r, sigma, M = 20, N = 20):
    dtau = T/N 

    S_partition = np.linspace(0.0, S_max, M+1)
    tau_partition = np.linspace(0.0, T, N+1)

    V = np.zeros((N+1, M+1))                     #V[n,i] = option value at time to maturity tau_n := tau_partition[n] and stock price S_i := S_partition[i]
    V[0, :] = np.maximum(K - S_partition, 0)            #initial condition at n = 0 (maturity)
    V[:, 0] = K * np.exp(-r * tau_partition)         #first boundary condition: if price of asset is zero, value of option is the discounted strike price
    V[:, -1] = 0.0                               #second boundary condition: we assume S_max is sufficiently large so that the option is essentially worthless if the asset reaches this price

    i = np.arange(1,M)      #vector of interior 'spatial' indices

    a = -0.5 * dtau * ((sigma**2 * i**2) - (r * i))
    b = 1.0 + dtau * ((sigma**2 * i**2) + r)
    c = -0.5 * dtau * ((sigma**2 * i**2) + (r * i))
    A = np.diag(b) + np.diag(c[:-1], 1) + np.diag(a[1:], -1)    #tridiagonal matrix for the finite difference scheme

    for n in range(N):
        rhs = V[n, 1:M].copy()
        rhs[0] -= a[0] * V[n+1, 0]
        rhs[-1] -= c[-1] * V[n+1, M]
        
        V[n+1, 1:M] = np.linalg.solve(A, rhs)
    
    return S_partition, tau_partition, V 


def european_put_closed_form(S, K, T, r, sigma):
    if T == 0: 
        return max(K - S, 0.0)
    if K == 0:
        return 0.0
    if S == 0:
        return K * np.exp(-r * T)
    
    N = NormalDist().cdf

    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)

    return K * np.exp(-r * T) * N(-d2) - S * N(-d1)


def european_put_closed_form_sampled(S_max, K, T, r, sigma, M = 20, N = 20):
    S_partition = np.linspace(0.0, S_max, M+1)
    tau_partition = np.linspace(0.0, T, N+1)

    V = np.zeros((N+1, M+1))      
    for i in range(N+1):
        for j in range(M+1):
            V[i,j] = european_put_closed_form(S_partition[j], K, tau_partition[i], r, sigma)    

    return S_partition, tau_partition, V      


def plot_option_surface(S_partition1, tau_partition1, V1, S_partition2 = None, tau_partition2 = None, V2 = None, stride = 1):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    S1, Tau1 = np.meshgrid(S_partition1, tau_partition1)
    ax.plot_wireframe(Tau1[::stride, ::stride], S1[::stride, ::stride], V1[::stride, ::stride], color='red', linewidth=0.6)
    red_patch = mpatches.Patch(color='red', label='Finite difference approximation')
    
    if V2 is not None:
        if S_partition2 is not None and tau_partition2 is not None:
            S2, Tau2 = np.meshgrid(S_partition2, tau_partition2)
            ax.plot_wireframe(Tau2[::stride, ::stride], S2[::stride, ::stride], V2[::stride, ::stride], color='blue', linewidth=0.6)
        else:
            ax.plot_wireframe(Tau1[::stride, ::stride], S1[::stride, ::stride], V2[::stride, ::stride], color='blue', linewidth=0.6)

        blue_patch = mpatches.Patch(color='blue', label='Exact solution (sampled on mesh vertices)')
        plt.legend(handles=[red_patch, blue_patch], loc = 'best')

    ax.set_xlabel("Time to maturity $\\tau$")
    ax.set_ylabel("Stock price $S$")
    ax.set_zlabel("European put option value $V$")
    ax.set_title("European put option value")
    plt.tight_layout()
    plt.show()


def plot_option_curve(S_partition1, V1, S_partition2 = None, V2 = None):
    fig = plt.figure(figsize=(10,6))

    plt.plot(S_partition1, V1[-1, :], label="Finite difference approximation", color = 'red')

    if V2 is not None:
        if S_partition2 is not None:
            plt.plot(S_partition2, V2[-1, :], label = "Exact solution (sampled on partition points)", color = 'blue')
        else:
            plt.plot(S_partition1, V2[-1, :], label = "Exact solution (sampled on partition points)")
    plt.xlabel("Stock price $S$")
    plt.ylabel("European put option value at $t=0$")
    plt.title("European put option value at $t = 0$")
    plt.legend(loc = 'best')
    plt.tight_layout()
    plt.show()

S_partition1, tau_partition1, V_fdm = european_put_implicit_fdm(400, 100, 1.0, 0.05, 0.2, 10, 10)
S_partition2, tau_partition2, V_exact = european_put_closed_form_sampled(400, 100, 1.0, 0.05, 0.2, 10, 10)

plot_option_surface(S_partition1, tau_partition1, V_fdm, S_partition2, tau_partition2, V_exact)
plot_option_curve(S_partition1, V_fdm, S_partition2, V_exact)




