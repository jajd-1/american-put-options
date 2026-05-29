import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import european_fdm as efdm
from pathlib import Path

script_dir = Path(__file__).resolve().parent.parent
images_dir = script_dir/'images'


def plot_option_surface(S_partition1, tau_partition1, V1, S_partition2 = None, tau_partition2 = None, V2 = None, stride = 1, 
                        label1 = '', label2 = '', xlabel = '', ylabel = '', zlabel = '', title = '', 
                        american = False, boundary_S = None, boundary_V = None, projection_height = 0, connector_step = 20):
    
    fig = plt.figure(figsize = (10, 8))
    ax = fig.add_subplot(111, projection = '3d')

    S1, Tau1 = np.meshgrid(S_partition1, tau_partition1)
    ax.plot_wireframe(Tau1, S1, V1, color='red', linewidth=0.6)

    if V2 is not None:
        if S_partition2 is not None and tau_partition2 is not None:
            S2, Tau2 = np.meshgrid(S_partition2, tau_partition2)
            ax.plot_wireframe(Tau2, S2, V2, color = 'blue', linewidth = 0.6)
        else:
            ax.plot_wireframe(Tau1, S1, V2, color = 'blue', linewidth = 0.6)

        if (label1 != '' and label2 != ''):
            red_patch = mpatches.Patch(color = 'red', label = label1)
            blue_patch = mpatches.Patch(color='blue', label = label2)
            plt.legend(handles = [red_patch, blue_patch], loc = 'best')
    
    if american == True: 
        if ((boundary_S is None) or (boundary_V is None)):
            raise ValueError
        
        valid = ~np.isnan(boundary_S) & ~np.isnan(boundary_V)
        tau_valid = tau_partition1[valid]
        boundary_S_valid = boundary_S[valid]
        boundary_V_valid = boundary_V[valid]

        ax.plot(tau_valid, boundary_S_valid, np.full_like(boundary_S_valid, projection_height), color = 'black', linestyle = '--', linewidth = 2, 
                label = 'Free boundary in $(S,\\tau)$ plane')

        ax.plot(tau_valid, boundary_S_valid, boundary_V_valid, color = 'blue', linewidth = 2.5, label = 'Free boundary projected onto surface')

        for tau_i, s_i, v_i in zip(tau_valid[::connector_step], boundary_S_valid[::connector_step], boundary_V_valid[::connector_step]):
            ax.plot([tau_i, tau_i], [s_i, s_i], [projection_height, v_i], color = 'gray', linewidth = 0.8, alpha = 0.7)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_title(title)
    plt.tight_layout()

    filename = title.replace(' ', '_') + '.png'
    plt.savefig(images_dir/filename, dpi = 300, bbox_inches = 'tight')
    
    plt.show()


def plot_option_curve(S_partition1, V1, S_partition2 = None, V2 = None, label1 = '', label2 = '', xlabel = '', ylabel = '', title = ''):
    fig = plt.figure(figsize = (10,6))

    plt.plot(S_partition1, V1[-1, :], label = label1, color = 'red')

    if V2 is not None:
        if S_partition2 is not None:
            plt.plot(S_partition2, V2[-1, :], label = label2, color = 'blue')
        else:
            plt.plot(S_partition1, V2[-1, :], label = label2, color = 'blue')
            
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    if (label1 != '' and label2 != ''):
        plt.legend(loc = 'best')

    plt.tight_layout()

    filename = title.replace(' ', '_') + '.png'
    plt.savefig(images_dir/filename, dpi = 300, bbox_inches = 'tight')

    plt.show()


def plot_errors_by_meshsize(start, stop, step):
    errors = {}

    for n in range(start, stop, step):
        _, _, V_fdm = efdm.european_put_implicit_fdm(400, 100, 1.0, 0.05, 0.2, n, n)             #change last entry to n*n for better mesh
        _, _, V_exact = efdm.european_put_closed_form_sampled(400, 100, 1.0, 0.05, 0.2, n, n)    #ditto
        errors[n] = np.max(np.abs(V_fdm - V_exact))

    fig = plt.figure(figsize = (10,6))
    plt.plot(errors.keys(), errors.values(), marker = '')
    plt.xlabel('Inverse mesh size')
    plt.ylabel('$L^\\infty$ error')
    plt.title('$L^\\infty$ error of finite difference approximation in terms of mesh size')
    plt.grid(True)

    filename = 'meshsize_errors.png'
    plt.savefig(images_dir/filename, dpi = 300, bbox_inches = 'tight')

    plt.show()


def plot_free_boundary(tau_partition, boundary_S):
    fig = plt.figure(figsize = (10, 6))
    plt.plot(tau_partition, boundary_S, color = 'black', linewidth = 2)
    plt.xlabel('Time to maturity $\\tau$')
    plt.ylabel('Exercise threshold')
    plt.title('Approximate free boundary for an American put option')
    plt.grid(True)
    plt.tight_layout()

    filename = 'approx_free_boundary.png'
    plt.savefig(images_dir/filename, dpi = 300, bbox_inches = 'tight')

    plt.show()


# ------ Family plots -------


def plot_option_curve_family(curve_dict, parameter_symbol = r'\sigma'):
    fig = plt.figure(figsize = (10,6))

    for param_value, (S_partition, _, V, _) in curve_dict.items():
        plt.plot(S_partition, V[-1, :], label = fr'${parameter_symbol} = {param_value}$')
            
    plt.xlabel('Underlying price $S$')
    plt.ylabel('Option value at $t=0$')
    plt.title(fr'Option value at $t=0$ for different ${parameter_symbol}$')
    plt.legend(loc = 'best')
    plt.tight_layout()

    filename = 'option_value_varying_' + parameter_symbol.replace('\\', '') + '.png'
    plt.savefig(images_dir/filename, dpi = 300, bbox_inches = 'tight')

    plt.show()


def plot_free_boundary_family(curve_dict, parameter_symbol = r'\sigma'):
    fig = plt.figure(figsize = (10, 6))

    for param_value, (_, tau_partition, _, boundary) in curve_dict.items():
        plt.plot(tau_partition, boundary, linewidth=2, label = fr'${parameter_symbol} = {param_value}$')

    plt.xlabel('Time to maturity $\\tau$')
    plt.ylabel('Exercise threshold')
    plt.title(fr'Approximate free boundary for an American put option for different ${parameter_symbol}$')
    plt.grid(True)
    plt.legend(loc = 'best')
    plt.tight_layout()

    filename = 'free_boundary_varying_' + parameter_symbol.replace('\\', '') + '.png'
    plt.savefig(images_dir/filename, dpi = 300, bbox_inches = 'tight')

    plt.show()
