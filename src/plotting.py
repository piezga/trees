# src/plotting.py
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import cmocean
import matplotlib.colors as mcolors
from matplotlib import patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from src.compute import marchenko_pastur_pdf, marchenko_pastur_bounds

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
    """
    plt.rcParams.update({
        'font.size': 14,
        'font.weight': 'bold',
        'axes.labelweight': 'bold',
        'axes.titlesize': 16,
        'axes.titleweight': 'bold',
        'lines.linewidth': 2.5,
        'lines.markersize': 8
    })
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))
    
    x = np.arange(1, num_species + 1)
    colors = sns.color_palette("muted", n_colors=2)
    
    # --- Main Plot (50 m) ---
    ax.plot(x, res_50_senm, 'o--', color=colors[0], label='SENM (50 m)')
    ax.fill_between(
        x, res_50_senm - res_50_senm_std, res_50_senm + res_50_senm_std,
        color=colors[0], alpha=0.25
    )
    for spectrum in res_50_forest_list:
        ax.plot(x, spectrum, 'o-', color=colors[0], alpha=0.5, markerfacecolor='white')
    
    # --- Inset (5 m) ---
    axins = inset_axes(ax, width="45%", height="45%", loc='lower left', borderpad=2)
    axins.plot(x, res_5_senm, 'o--', color=colors[1], label='SENM (5 m)')
    axins.fill_between(
        x, res_5_senm - res_5_senm_std, res_5_senm + res_5_senm_std,
        color=colors[1], alpha=0.25
    )
    for spectrum in res_5_forest_list:
        axins.plot(x, spectrum, 'o-', color=colors[1], alpha=0.5, markerfacecolor='white')
    
    # --- Axis Formatting ---
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1, 90)
    ax.set_ylim(1e-3, 100)
    ax.set_xlabel(r'$\lambda$ Rank', labelpad=8)
    ax.set_ylabel(r'$\lambda$ Value', labelpad=8)
    ax.set_title('Spectral Comparison at 50 m ')
    ax.grid(True, which="both", linestyle='--', alpha=0.5)
    
    axins.set_xscale('log')
    axins.set_yscale('log')
    axins.grid(True, linestyle='--', alpha=0.5)
    axins.set_xlim(1, num_species)
    axins.set_ylim(0.5, max(np.max(res_5_senm + res_5_senm_std),
                            max(np.max(spectrum) for spectrum in res_5_forest_list)))
    axins.set_title('5 m ', fontsize=12, pad=2)
    axins.tick_params(labelsize=8)
    
    # --- Label and Layout ---
    ax.text(-0.12, 1.05, '(a)', transform=ax.transAxes, fontsize=24, weight='bold')
    ax.legend(loc='best', frameon=False, fontsize=12)
    plt.tight_layout()
    
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
    """
    plt.rcParams.update({
        'font.size': 14,
        'font.weight': 'bold',
        'axes.labelweight': 'bold',
        'axes.titlesize': 16,
        'axes.titleweight': 'bold',
        'lines.linewidth': 2.5,
        'lines.markersize': 8
    })
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))
    
    colors_b = sns.color_palette("colorblind", n_colors=len(species_array))
    x = np.array(resolutions_list, dtype=float)
    
    for idx, num_species in enumerate(species_array):
        ax.plot(
            x,
            comm_diffs_dict[num_species],
            'o--',
            label=f"{num_species} Species",
            color=colors_b[idx],
            markerfacecolor='white',
            linewidth=2
        )
    
    # --- Axis labels + styling ---
    ax.set_xlabel('Inverse resolution (m)', labelpad=8)
    ax.set_ylabel(r'$\Delta N_c$', labelpad=8)
    ax.set_title('Difference in detected communities')
    
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(fontsize=12, loc='best', frameon=False)
    ax.text(-0.12, 1.05, '(b)', transform=ax.transAxes, fontsize=24, weight='bold')
    plt.tight_layout()
    
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
    
    plt.xlabel('Inverse resolution (m)')
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

def plot_correlation_matrices_comparison(
    unfiltered_corr,
    filtered_reordered_corr,
    data_type='SENM',
    cmap=None,
    vmin=-0.5,
    vmax=0.5,
    figsize=(10, 4),
    filename=None,
    show=True
):
    """
    Plot side-by-side comparison of unfiltered and filtered/reordered correlation matrices.
    
    Parameters
    ----------
    unfiltered_corr : array
        Unfiltered correlation matrix
    filtered_reordered_corr : array
        Filtered and community-reordered correlation matrix
    data_type : str
        Label for the data type (e.g., 'SENM', 'Forest')
    cmap : matplotlib colormap, optional
        Colormap to use. Defaults to cmocean.cm.balance
    vmin, vmax : float
        Color scale limits
    figsize : tuple
        Figure size (width, height)
    filename : str, optional
        Path to save figure
    show : bool
        Whether to display the plot
        
    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : array of matplotlib.axes.Axes
    """
    if cmap is None:
        cmap = cmocean.cm.balance
    
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    
    fig, axes = plt.subplots(1, 2, figsize=figsize, constrained_layout=True)
    
    # Unfiltered correlation matrix
    axes[0].imshow(unfiltered_corr, cmap=cmap, norm=norm)
    axes[0].set_title(f"{data_type} – Unfiltered")
    axes[0].set_xlabel("Species index")
    axes[0].set_ylabel("Species index")
    
    # Filtered & reordered correlation matrix
    im1 = axes[1].imshow(filtered_reordered_corr, cmap=cmap, norm=norm)
    axes[1].set_title(f"{data_type} – Filtered & Reordered by Community")
    axes[1].set_xlabel("Species index")
    axes[1].set_ylabel("Species index")
    
    # Colorbar
    fig.colorbar(im1, ax=axes, orientation="vertical", fraction=0.03, 
                 pad=0.04, label="Correlation")
    fig.suptitle(f"{data_type} Correlation Matrices", fontsize=14, weight="bold")
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig, axes

from scipy.cluster.hierarchy import dendrogram

def plot_community_dendrogram(
    linkage_matrix,
    threshold=1e-4,
    ax=None,
    orientation='top',
    leaf_rotation=0,
    leaf_font_size=10,
    color_threshold_color='#ED2939',
    figsize=(8, 6),
    normalize=True,
    title=None,
    ylabel=r'$D/D_{max}$',
    xlabel='Node index',
    ylim_padding=(0.2, 0.2),  # Changed to symmetric padding
    ylim_auto=True,  # New parameter
    ylim_manual=None,  # New parameter for manual override
    filename=None,
    show=True
):
    """
    Plot hierarchical clustering dendrogram with custom styling.
    
    Parameters
    ----------
    linkage_matrix : array
        Linkage matrix from scipy.cluster.hierarchy.linkage
    threshold : float
        Distance threshold for cutting the dendrogram (horizontal line)
    ax : matplotlib.axes.Axes, optional
        Axes to plot on
    orientation : str
        Dendrogram orientation ('top', 'bottom', 'left', 'right')
    leaf_rotation : float
        Rotation angle for leaf labels
    leaf_font_size : int
        Font size for leaf labels
    color_threshold_color : str
        Color for the threshold line
    figsize : tuple
        Figure size if ax is None
    normalize : bool
        If True, normalize linkage distances by maximum distance
    title : str, optional
        Plot title
    ylabel : str
        Y-axis label
    xlabel : str
        X-axis label
    ylim_padding : tuple
        (bottom_padding_fraction, top_padding_fraction) for y-axis limits
    ylim_auto : bool
        If True, automatically determine y-limits based on data range
    ylim_manual : tuple, optional
        Manual (ymin, ymax) override. Takes precedence over ylim_auto
    filename : str, optional
        Path to save figure
    show : bool
        Whether to display the plot
        
    Returns
    -------
    ax : matplotlib.axes.Axes
    dendrogram_dict : dict
        Dictionary of data structures computed to render the dendrogram
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    
    # Normalize linkage matrix if requested
    if normalize:
        linkage_normalized = linkage_matrix.copy()
        tmax = linkage_matrix[:, 2][-1]
        linkage_normalized[:, 2] = linkage_matrix[:, 2] / tmax
    else:
        linkage_normalized = linkage_matrix
    
    # Create labels
    n_samples = linkage_matrix.shape[0] + 1
    labelList = [i + 1 for i in range(n_samples)]
    
    # Plot dendrogram
    dendrogram_dict = dendrogram(
        linkage_normalized,
        labels=labelList,
        ax=ax,
        leaf_rotation=leaf_rotation,
        orientation=orientation,
        color_threshold=threshold,
        above_threshold_color='k',
        leaf_font_size=leaf_font_size
    )
    
    # Set y-axis limits
    if ylim_manual is not None:
        # Manual override
        ax.set_ylim(ylim_manual)
    elif ylim_auto:
        # Automatic: base on actual data range
        distances = linkage_normalized[:, 2]
        dmin = distances[0]  # Minimum distance
        dmax = distances[-1]  # Maximum distance
        
        # Use log-scale aware padding
        # In log space, we want to add/subtract a fraction of the log range
        log_dmin = np.log10(max(dmin, 1e-10))  # Avoid log(0)
        log_dmax = np.log10(dmax)
        log_range = log_dmax - log_dmin
        
        # Add padding in log space
        ymin = 10 ** (log_dmin - ylim_padding[0] * log_range)
        ymax = 10 ** (log_dmax + ylim_padding[1] * log_range)
        
        ax.set_ylim(ymin, ymax)
    else:
        # Old behavior: use fixed padding
        tmin = linkage_normalized[:, 2][0] - ylim_padding[0] * linkage_normalized[:, 2][0]
        tmax = linkage_normalized[:, 2][-1] + ylim_padding[1] * linkage_normalized[:, 2][-1]
        ax.set_ylim(tmin, tmax)
    
    # Add threshold line
    ax.axhline(y=threshold, color=color_threshold_color, linestyle='--', linewidth=2, 
               label=f'Cut threshold')
    
    # Labels and styling
    if title:
        ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_yscale('log')
    
    # Set y-ticks dynamically based on the actual range
    current_ylim = ax.get_ylim()
    # Generate reasonable tick positions in log scale
    log_ymin = np.floor(np.log10(current_ylim[0]))
    log_ymax = np.ceil(np.log10(current_ylim[1]))
    yticks = [10**i for i in range(int(log_ymin), int(log_ymax) + 1)]
    ax.set_yticks(yticks)
    
    ax.set_xticks([])  # Remove x-ticks as in the example
    ax.legend(loc='best', fontsize=10)
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    
    return ax, dendrogram_dict

def supplementary_plots(results_for_supplementary, verbose):
    # Module that generates supplementary_plots
    # Needs the results from the first subplot
    for result in results_for_supplementary:
        num_species = result['num_species']
        resolution = result['resolution']

        # Plot spectra with MP bounds
        plot_spectra_with_mp_bounds(
            result['forest_spectrum'],
            result['senm_spectrum'],
            result['lambda_max_forest'],
            result['lambda_max_senm'],
            resolution,
            num_species,
            filename=f'spectra/{num_species}_species_resolution_{resolution}.png',
            show=False,
            verbose=verbose
        )

        # Plot eigenvalue density
        plot_eigenvalue_density_vs_mp(
            result['forest_spectrum'],
            result['senm_spectrum'],
            result['bins'],
            resolution,
            num_species,
            result['lambda_max_forest'],
            result['lambda_max_senm'],
            filename=f'spectra/{num_species}_species_density_resolution_{resolution}.png',
            show=False,
            verbose=verbose
        )

    # Plot community count vs resolution for each species count
    # Group results by num_species
    from collections import defaultdict
    community_data = defaultdict(lambda: {'resolutions': [], 'forest': [], 'senm': []})

    for result in results_for_supplementary:
        ns = result['num_species']
        community_data[ns]['resolutions'].append(result['resolution'])
        community_data[ns]['forest'].append(result['forest_communities'])
        community_data[ns]['senm'].append(result['senm_communities'])

    for num_species, data in community_data.items():
        plot_community_count_vs_resolution(
            data['resolutions'],
            data['forest'],
            data['senm'],
            num_species,
            filename=f'spectra/{num_species}_communities_vs_resolution.png',
            show=False,
            verbose=verbose
        )

def plot_community_comparison(
    corr_matrix,
    CM_method1,
    CM_method2,
    method1_name='Laplacian Diffusion',
    method2_name='Ward Clustering',
    cmap=None,
    vmin=-0.5,
    vmax=0.5,
    figsize=(14, 6),
    filename=None,
    show=True
):
    """
    Compare two community detection methods side-by-side.
    
    Parameters
    ----------
    corr_matrix : array
        Original correlation matrix
    CM_method1 : array
        Community membership from first method
    CM_method2 : array
        Community membership from second method
    method1_name : str
        Name of first method
    method2_name : str
        Name of second method
    cmap : matplotlib colormap, optional
        Colormap for correlation matrix
    vmin, vmax : float
        Color scale limits
    figsize : tuple
        Figure size
    filename : str, optional
        Path to save figure
    show : bool
        Whether to display the plot
        
    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : array of matplotlib.axes.Axes
    """
    from src.compute import get_community_bounds

    if cmap is None:
        cmap = cmocean.cm.balance
    
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    
    fig, axes = plt.subplots(1, 2, figsize=figsize, constrained_layout=True)
    
    # Method 1: Reorder by communities
    idx1 = np.argsort(CM_method1)
    reordered1 = np.array([[corr_matrix[i][j] for j in idx1] for i in idx1])
    
    im1 = axes[0].imshow(reordered1, cmap=cmap, norm=norm)
    axes[0].set_title(f"{method1_name}\n({len(np.unique(CM_method1))} communities)")
    axes[0].set_xlabel("Species index")
    axes[0].set_ylabel("Species index")
    
    # Add community boundaries for method 1
    bounds1 = get_community_bounds(CM_method1[idx1])
    draw_community_boxes(axes[0], bounds1)
    
    # Method 2: Reorder by communities
    idx2 = np.argsort(CM_method2)
    reordered2 = np.array([[corr_matrix[i][j] for j in idx2] for i in idx2])
    
    im2 = axes[1].imshow(reordered2, cmap=cmap, norm=norm)
    axes[1].set_title(f"{method2_name}\n({len(np.unique(CM_method2))} communities)")
    axes[1].set_xlabel("Species index")
    axes[1].set_ylabel("Species index")
    
    # Add community boundaries for method 2
    bounds2 = get_community_bounds(CM_method2[idx2])
    draw_community_boxes(axes[1], bounds2)
    
    # Shared colorbar
    fig.colorbar(im2, ax=axes, orientation="vertical", 
                 fraction=0.03, pad=0.04, label="Correlation")
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig, axes


def draw_community_boxes(ax, bounds):
    """Draw boxes around communities on a matrix plot."""
    bounds = np.array(bounds)
    bounds[0] += 0.2  # Small offset for visibility
    bounds[-1] -= 0.2
    
    for n, edge in enumerate(np.diff(bounds)):
        ax.add_patch(patches.Rectangle(
            (bounds[n], bounds[n]),
            edge, edge, 
            fill=False, 
            linewidth=1.5, 
            ls='--',
            edgecolor='black'
        ))




def plot_community_confusion_matrix(
    CM_method1,
    CM_method2,
    method1_name='Method 1',
    method2_name='Method 2',
    figsize=(8, 6),
    filename=None,
    show=True
):
    """
    Plot confusion matrix showing overlap between two community assignments.
    
    Parameters
    ----------
    CM_method1 : array
        Community membership from first method
    CM_method2 : array
        Community membership from second method
    method1_name : str
        Name of first method
    method2_name : str
        Name of second method
    figsize : tuple
        Figure size
    filename : str, optional
        Path to save figure
    show : bool
        Whether to display the plot
        
    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    ari : float
        Adjusted Rand Index
    """
    from src.compute import compute_community_overlap

    ari, confusion = compute_community_overlap(CM_method1, CM_method2)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    im = ax.imshow(confusion, cmap='Blues', aspect='auto')
    
    # Add text annotations
    for i in range(confusion.shape[0]):
        for j in range(confusion.shape[1]):
            text = ax.text(j, i, confusion[i, j],
                          ha="center", va="center", color="black")
    
    ax.set_xlabel(f"{method2_name} Communities")
    ax.set_ylabel(f"{method1_name} Communities")
    ax.set_title(f"Community Overlap (ARI = {ari:.3f})")
    
    # Set ticks
    ax.set_xticks(np.arange(confusion.shape[1]))
    ax.set_yticks(np.arange(confusion.shape[0]))
    ax.set_xticklabels(np.arange(1, confusion.shape[1] + 1))
    ax.set_yticklabels(np.arange(1, confusion.shape[0] + 1))
    
    plt.colorbar(im, ax=ax, label="Number of species")
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig, ax, ari
