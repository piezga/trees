# src/plotting.py
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from src.functions import marchenko_pastur_pdf, marchenko_pastur_bounds

def plot_dual_spectra_with_inset(
    res_50_senm, res_50_senm_std, res_50_forest_list,
    res_5_senm, res_5_senm_std, res_5_forest_list,
    num_species=100,
    ax=None,
    show=True,
    filename=None
):
    """
    Plot dual-resolution spectra comparison with inset.
    
    Parameters
    ----------
    res_50_senm : array
        SENM spectrum at 50m resolution
    res_50_senm_std : array
        Standard deviation for SENM at 50m
    res_50_forest_list : list of arrays
        Forest spectra at 50m resolution
    res_5_senm : array
        SENM spectrum at 5m resolution
    res_5_senm_std : array
        Standard deviation for SENM at 5m
    res_5_forest_list : list of arrays
        Forest spectra at 5m resolution
    num_species : int
        Number of species
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure
    show : bool
        Whether to call plt.show()
    filename : str, optional
        If provided, save figure to this path
        
    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes object with the plot
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))
    
    x = np.arange(1, num_species + 1)
    colors = sns.color_palette("muted", n_colors=2)
    
    # Main plot: 50m resolution
    ax.plot(x, res_50_senm, 'o--', color=colors[0], label='SENM (50 m)')
    ax.set_title('50 m', fontsize=22, pad=4)
    ax.fill_between(x, res_50_senm - res_50_senm_std, 
                     res_50_senm + res_50_senm_std,
                     color=colors[0], alpha=0.25)
    
    for spectrum in res_50_forest_list:
        ax.plot(x, spectrum, 'o-', color=colors[0], 
                alpha=0.5, markerfacecolor='white')
    
    # Inset: 5m resolution
    axins = inset_axes(ax, width="45%", height="45%", 
                       loc='lower left', borderpad=2)
    axins.plot(x, res_5_senm, 'o--', color=colors[1], label='SENM (5 m)')
    axins.fill_between(x, res_5_senm - res_5_senm_std, 
                       res_5_senm + res_5_senm_std,
                       color=colors[1], alpha=0.25)
    
    for spectrum in res_5_forest_list:
        axins.plot(x, spectrum, 'o-', color=colors[1], 
                   alpha=0.5, markerfacecolor='white')
    
    # Formatting
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1, 90)
    ax.set_ylim(1e-3, 100)
    ax.set_xlabel(r'$\lambda$ rank')
    ax.set_ylabel(r'$\lambda$')
    ax.grid(True, which="both", linestyle='--', alpha=0.5)
    
    axins.grid(True, linestyle='--', alpha=0.5)
    axins.set_xscale('log')
    axins.set_yscale('log')
    axins.set_xlim(1, num_species)
    axins.set_ylim(0.5, max(np.max(res_5_senm + res_5_senm_std),
                            max(np.max(spectrum) for spectrum in res_5_forest_list)))
    axins.set_title('5 m', fontsize=18, pad=4)
    axins.tick_params(labelsize=8)
    
    ax.text(-0.12, 1.05, '(a)', transform=ax.transAxes, 
            fontsize=24, weight='bold')
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    
    return ax

def plot_size_effect_panel(
    resolutions_list,
    comm_diffs_dict,
    species_array,
    ax=None,
    show=False
):
    """
    Plot Panel B: Community difference vs resolution for multiple species counts.
    
    Parameters
    ----------
    resolutions_list : list
        List of resolution values
    comm_diffs_dict : dict
        Dictionary mapping num_species -> list of community differences
    species_array : list
        List of species counts to plot
    ax : matplotlib.axes.Axes, optional
        Axes to plot on
    show : bool
        Whether to display the plot
        
    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))
    
    colors_b = sns.color_palette("colorblind", n_colors=len(species_array))
    x = np.array(resolutions_list)
    
    for idx, num_species in enumerate(species_array):
        ax.plot(
            x,
            comm_diffs_dict[num_species],
            'o--',
            label=f"N = {num_species}",
            color=colors_b[idx],
            markerfacecolor='white',
            linewidth=2
        )
    
    ax.set_xlabel("Resolution (m)")
    ax.set_ylabel(r"$\Delta N_c$")
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(fontsize=22, loc='best', frameon=False)
    ax.text(-0.12, 1.05, '(b)', transform=ax.transAxes, fontsize=24, weight='bold')
    
    if show:
        plt.show()
    
    return ax

def plot_spectra_with_mp_bounds(
    forest_spectrum,
    senm_spectrum,
    lambda_max_forest,
    lambda_max_senm,
    resolution,
    num_species,
    filename=None,
    show=False,
    verbose=False
):
    """
    Plot forest and SENM spectra with Marchenko-Pastur bounds.
    
    Parameters
    ----------
    forest_spectrum : array
        Forest eigenvalue spectrum
    senm_spectrum : array
        SENM eigenvalue spectrum
    lambda_max_forest : float
        MP maximum bound for forest
    lambda_max_senm : float
        MP maximum bound for SENM
    resolution : int
        Spatial resolution in meters
    num_species : int
        Number of species
    filename : str, optional
        Path to save figure
    show : bool
        Whether to display the plot
    verbose : bool
        Print save confirmation
        
    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    fig = plt.figure(figsize=(12, 7))
    
    plt.rcParams.update({
        'font.size': 14,
        'font.weight': 'bold',
        'axes.labelweight': 'bold',
        'axes.titlesize': 16,
        'axes.titleweight': 'bold',
        'lines.linewidth': 2.5,
        'lines.markersize': 8
    })
    
    x = np.arange(1, num_species + 1)
    
    plt.loglog(x, forest_spectrum, 'o-', label='Forest')
    plt.loglog(x, senm_spectrum, 'o-', label='SENM')
    
    plt.title(f'{resolution} m')
    plt.grid(True, alpha=0.5)
    plt.axhline(lambda_max_forest, label='MP Forest')
    plt.axhline(lambda_max_senm, color='orange', label='MP SENM')
    plt.legend()
    plt.ylim(1e-2, 1e2)
    
    if filename:
        plt.savefig(filename)
        if verbose:
            print(f'Saved spectra to {filename}')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig


def plot_eigenvalue_density_vs_mp(
    forest_spectrum,
    senm_spectrum,
    bins,
    resolution,
    num_species,
    lambda_max_forest,
    lambda_max_senm,
    filename=None,
    show=False,
    verbose=False
):
    """
    Plot eigenvalue density histograms with Marchenko-Pastur PDF overlay.
    
    Parameters
    ----------
    forest_spectrum : array
        Forest eigenvalue spectrum
    senm_spectrum : array
        SENM eigenvalue spectrum
    bins : array
        [n_bins_x_senm, n_bins_y_senm, n_bins_x_forest, n_bins_y_forest]
    resolution : int
        Spatial resolution in meters
    num_species : int
        Number of species
    lambda_max_forest : float
        MP maximum bound for forest
    lambda_max_senm : float
        MP maximum bound for SENM
    filename : str, optional
        Path to save figure
    show : bool
        Whether to display the plot
    verbose : bool
        Print save confirmation
        
    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    fig = plt.figure(figsize=(8, 6))
    
    # Histogram of eigenvalues
    plt.hist(forest_spectrum, bins=20, density=True, alpha=0.5,
             label='Forest eigenvalues', color='C0')
    plt.hist(senm_spectrum, bins=20, density=True, alpha=0.5,
             label='SENM eigenvalues', color='C1')
    
    # MP theoretical PDF
    q_forest = (bins[2] * bins[3]) / num_species
    q_senm = (bins[0] * bins[1]) / num_species
    l_vals = np.linspace(0, max(lambda_max_forest, lambda_max_senm) * 1.2, 400)
    
    mp_pdf_forest, lmin_forest, lmax_forest = marchenko_pastur_pdf(l_vals, q_forest)
    mp_pdf_senm, lmin_senm, lmax_senm = marchenko_pastur_pdf(l_vals, q_senm)
    
    plt.plot(l_vals, mp_pdf_forest, 'C0-', lw=2, label='MP Forest')
    plt.plot(l_vals, mp_pdf_senm, 'C1--', lw=2, label='MP SENM')
    
    plt.xlabel("Eigenvalue")
    plt.ylabel("Density")
    plt.title(f"Eigenvalue Spectrum Density vs MP ({resolution} m, N={num_species})")
    plt.legend()
    plt.grid(alpha=0.4)
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename)
        if verbose:
            print(f'Saved MP density plot to {filename}')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig


def plot_community_count_vs_resolution(
    resolutions,
    forest_communities,
    senm_communities,
    num_species,
    filename=None,
    show=False,
    verbose=False
):
    """
    Plot number of communities vs resolution for forest and SENM.
    
    Parameters
    ----------
    resolutions : array
        Resolution values in meters
    forest_communities : array
        Number of forest communities at each resolution
    senm_communities : array
        Number of SENM communities at each resolution
    num_species : int
        Number of species
    filename : str, optional
        Path to save figure
    show : bool
        Whether to display the plot
    verbose : bool
        Print save confirmation
        
    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    fig = plt.figure(figsize=(10, 6))
    
    plt.rcParams.update({
        'font.size': 14,
        'font.weight': 'bold',
        'axes.labelweight': 'bold',
        'axes.titlesize': 16,
        'axes.titleweight': 'bold',
        'lines.linewidth': 2.5,
        'lines.markersize': 8
    })
    
    x = np.array(resolutions, dtype=float)
    colors = sns.color_palette("colorblind", n_colors=2)
    
    plt.plot(x, forest_communities, 'o--', label='Forest Communities', color=colors[0])
    plt.plot(x, senm_communities, 'o--', label='SENM Communities', color=colors[1])
    
    plt.xlabel('Resolution (m)')
    plt.ylabel('Number of Communities')
    plt.title(f'Community Count vs Resolution ({num_species} species)')
    plt.grid(True, alpha=0.5)
    plt.legend()
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=300)
        if verbose:
            print(f'Saved community comparison plot to {filename}')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig


