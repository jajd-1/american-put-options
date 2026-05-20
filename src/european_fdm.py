import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches

def european_put_implicit_fdm(S_max, K, T, r, sigma, M = 500, N = 500):
    dS = S_max/M
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


def plot_option_surface(S_partition, tau_partition, V, stride = 1):
    S, Tau = np.meshgrid(S_partition, tau_partition)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_wireframe(Tau[::stride, ::stride], S[::stride, ::stride], V[::stride, ::stride], color='red', linewidth=0.6)

    ax.set_xlabel("Time to maturity $\\tau$")
    ax.set_ylabel("Stock price $S$")
    ax.set_zlabel("European put option value $V$")
    ax.set_title("European put option value computed using implicit finite differences")

    plt.tight_layout()
    plt.show()

def plot_option_curve(S_partition, V):
    plt.figure(figsize=(10,6))
    plt.plot(S_partition, V[-1, :], label="Finite difference")
    plt.xlabel("Stock price $S$")
    plt.ylabel("European put option value at $t=0$")
    plt.title("European put option value at $t = 0$ computed using implicit finite differences")
    plt.tight_layout()
    plt.show()

S_partition, tau_partition, V = european_put_implicit_fdm(400, 100, 1.0, 0.05, 0.2)
plot_option_surface(S_partition, tau_partition, V)
plot_option_curve(S_partition, V)




