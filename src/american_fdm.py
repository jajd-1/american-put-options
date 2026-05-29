import numpy as np


def american_put_implicit_psor(S_max, K, T, r, sigma, M = 500, N = 500, omega = 1.2, tol = 1e-8, max_iter = 10_000):
    dtau = T/N 

    S_partition = np.linspace(0, S_max, M+1)
    tau_partition = np.linspace(0, T, N+1)

    V = np.zeros((N+1, M+1))                            #V[n,i] = option value at time to maturity tau_n := tau_partition[n] and stock price S_i := S_partition[i]
    V[0, :] = np.maximum(K - S_partition, 0)            #initial condition at n = 0 (maturity)
    V[:, 0] = K                                         #first boundary condition: if price of asset is zero, we can exercise immediately at the strike price
    V[:, -1] = 0.                                       #second boundary condition: we assume S_max is sufficiently large so that the option is essentially worthless if the asset reaches S_max

    i = np.arange(1,M)      #vector of interior 'spatial' indices

    a = -0.5 * dtau * ((sigma**2 * i**2) - (r * i))
    b = 1 + dtau * ((sigma**2 * i**2) + r)
    c = -0.5 * dtau * ((sigma**2 * i**2) + (r * i))

    for n in range(N):
        rhs = V[n, 1:M].copy()
        rhs[0] -= a[0] * V[n+1, 0]
        rhs[-1] -= c[-1] * V[n+1, M]            #this does nothing as right boundary value is 0

        obstacle_vector = np.maximum(K - S_partition[1:M], 0) 

        x = V[n, 1:M].copy()        #initial guess for PSOR 

        for _ in range(max_iter):
            x_old = x.copy()
            
            for j in range(M-1):
                if j == 0:
                    left = V[n+1, 0]
                else:
                    left = x[j-1]
                
                if j == M-2:
                    right = V[n+1, M]
                else:
                    right = x_old[j+1]
                
                gauss_seidel_update = (rhs[j] - a[j]*left - c[j]*right) / b[j]

                over_relaxation_update = x_old[j] + omega*(gauss_seidel_update - x_old[j])

                x[j] = max(over_relaxation_update, obstacle_vector[j])
            
            if np.max(np.abs(x - x_old)) < tol:
                break

        V[n+1, 1:M] = x 
    
    return S_partition, tau_partition, V 


def extract_free_boundary(S_partition, V, K, tol=1e-6):
    barrier = np.maximum(K - S_partition, 0)

    boundary_S = np.full(V.shape[0], np.nan)
    boundary_V = np.full(V.shape[0], np.nan)

    for n in range(V.shape[0]):
        premium = V[n, :] - barrier

        exercise_indices = np.where((S_partition < K) & (premium <= tol))[0]

        if len(exercise_indices) > 0:
            i = exercise_indices[-1]
            boundary_S[n] = S_partition[i]
            boundary_V[n] = V[n, i]

    return boundary_S, boundary_V


