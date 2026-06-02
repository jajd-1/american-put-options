import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
from statistics import NormalDist
import european_fdm as efdm
import american_fdm as afdm 
import plotting_functions as plfn
import greeks as gr

S_partition1, tau_partition1, V_fdm = efdm.european_put_implicit_fdm(400, 100, 1.0, 0.05, 0.2, 1000, 1000)
S_partition2, tau_partition2, V_exact = efdm.european_put_closed_form_sampled(400, 100, 1.0, 0.05, 0.2, 1000, 1000)
S_partition3, tau_partition3, V_psor = afdm.american_put_implicit_psor(400, 100, 1.0, 0.05, 0.2, 1000, 1000)

boundary_S, boundary_V = afdm.extract_free_boundary(S_partition3, V_psor, 100)

Delta_am, Gamma_am = gr.compute_greeks(S_partition3, V_psor)
Delta_eur, Gamma_eur = gr.compute_greeks(S_partition1, V_fdm) 

boundary_Delta_am = np.full_like(boundary_S, np.nan, dtype = float)
boundary_Gamma_am = np.full_like(boundary_S, np.nan, dtype = float)

for n in range(len(boundary_S)):
    if not np.isnan(boundary_S[n]):
        i = np.where(S_partition3 == boundary_S[n])[0]
        if len(i) > 0:
            boundary_Delta_am[n] = Delta_am[n, i[0]]
            boundary_Gamma_am[n] = Gamma_am[n, i[0]]


plfn.plot_option_surface(S_partition1, tau_partition1, V_fdm, S_partition2, tau_partition2, V_exact,
                    label1 = 'Finite difference approximation',
                    label2 = 'Exact solution (sampled on mesh vertices)',
                    xlabel = 'Time to maturity $\\tau$',
                    ylabel = 'Underlying price $S$',
                    zlabel = 'European put option value $V$', 
                    title = 'European put option value as a function of the S and tau - FDM vs exact solution')


plfn.plot_option_curve(S_partition1, V_fdm, S_partition2, V_exact,
                    label1 = 'Finite difference approximation',
                    label2 = 'Exact solution (sampled on partition points)',
                    xlabel = 'Underlying price $S$',
                    ylabel = 'European put option value at $t=0$',
                    title = 'European put option value at $t=0$ as a function of S - FDM vs exact solution')


#plfn.plot_errors_by_meshsize(20, 800, 8)


plfn.plot_free_boundary(tau_partition3, boundary_S)

plfn.plot_near_expiry_asymptotic(tau_partition3, boundary_S, 100, 0.2, tau_min=0.0, tau_max=0.1)



plfn.plot_option_curve(S_partition1, V_fdm, S_partition3, V_psor,
                    label1 = 'European',
                    label2 = 'American',
                    xlabel = 'Underlying $S$',
                    ylabel = 'Option value at $t=0$',
                    title = 'American (PSOR) vs European (FDM) put option values at $t=0$')


plfn.plot_option_curve(S_partition1, V_psor - V_fdm,
                    xlabel = 'Underlying price $S$',
                    ylabel = 'Price difference at $t=0$ (American minus European)',
                    title = 'American (PSOR) minus European (FDM) put option value at $t=0$', 
                    vline = boundary_S[-1])


plfn.plot_option_surface(S_partition1, tau_partition1, V_psor, 
                    xlabel = 'Time to maturity $\\tau$',
                    ylabel = 'Underlying price $S$',
                    zlabel = 'American put option value $V$', 
                    title = 'American put option value as a function of S and tau, with free boundary', 
                    american = True, boundary_S = boundary_S, boundary_V = boundary_V)


plfn.plot_option_surface(S_partition1, tau_partition1, V_psor - V_fdm, view = (22, 146, 0),
                    xlabel = 'Time to maturity $\\tau$',
                    ylabel = 'Underlying price $S$',
                    zlabel = 'Price difference (American minus European)', 
                    title = 'American (PSOR) minus European (FDM) put option value as a function of S and tau')


plfn.plot_option_curve(S_partition1, Delta_eur, S_partition3, Delta_am,
                       label1 = 'European',
                       label2 = 'American',
                       xlabel = 'Underlying price $S$',
                       ylabel = 'Delta at $t=0$', 
                       title = 'Delta of European (FDM) and American (PSOR) put option at $t=0$', 
                       vline = boundary_S[-1])


plfn.plot_option_surface(S_partition3, tau_partition3, Delta_am, 
                        xlabel = 'Time to maturity $\\tau$',
                        ylabel = 'Underlying price $S$',
                        zlabel = 'Delta', 
                        title = 'Delta of American put option (PSOR) as function of $S$ and $tau$, with free boundary',
                        american = True, boundary_S = boundary_S, boundary_V = boundary_Delta_am)


plfn.plot_option_curve(S_partition3, Gamma_eur, S_partition1, Gamma_am,
                       label1 = 'European',
                       label2 = 'American',
                       xlabel = 'Underlying price $S$',
                       ylabel = 'Gamma at $t=0$', 
                       title = 'Gamma of European (FDM) and American (PSOR) put option at $t=0$', 
                       vline = boundary_S[-1])


plfn.plot_option_surface(S_partition3, tau_partition3[200:], Gamma_am[200:, :], 
                        xlabel = 'Time to maturity $\\tau$',
                        ylabel = 'Underlying price $S$',
                        zlabel = 'Gamma', 
                        title = 'Gamma of American put option (PSOR) as function of $S$ and $tau$')


# -------- Parameter variation ---------

sigmas = [0.1, 0.2, 0.3, 0.4]
boundary_dict = {}

for sigma in sigmas:
    S_partition, tau_partition, V = afdm.american_put_implicit_psor(400, 100, 1.0, 0.05, sigma, 500, 500)
    boundary_S, _ = afdm.extract_free_boundary(S_partition, V, 100)
    boundary_dict[sigma] = (S_partition, tau_partition, V, boundary_S)

plfn.plot_free_boundary_family(boundary_dict, parameter_symbol = r'\sigma')
plfn.plot_option_curve_family(boundary_dict, parameter_symbol = r'\sigma')


rates = [0.01, 0.03, 0.05, 0.07]
boundary_dict = {}

for r in rates:
    S_partition, tau_partition, V = afdm.american_put_implicit_psor(400, 100, 1.0, r, 0.2, 500, 500)
    boundary_S, _ = afdm.extract_free_boundary(S_partition, V, 100)
    boundary_dict[r] = (S_partition, tau_partition, V, boundary_S)

plfn.plot_free_boundary_family(boundary_dict, parameter_symbol = 'r')
plfn.plot_option_curve_family(boundary_dict, parameter_symbol = 'r')
