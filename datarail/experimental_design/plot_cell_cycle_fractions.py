import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import itertools
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def plot_stacked_bar(dfi, plot_num=None,
                     data_cols=['G1', 'S', 'G2', 'S_dropout',
                                'subG1', 'beyondG2', 'M', 'increase_fraction_dead'],
                     title=None,
                     show_legend=False):
    """Make stacked barplot for a given cell line and drug.
    
    Parameters
    ----------
    data : pandas.DataFrame 
          report dataframe from the 'run_cell_cycle_gating' function
    plot_num : Optional[int]
       Plot number in the context of mulitple such plots.
       Default is None.
    data_cols : Optional[list of str]
       Columns to be used to make stacked barplots
       Default it is ['G1', 'S', 'G2', 'S_dropout', 'other', 'M', 'fraction_dead'].
    title : Optional[str]
       Name or title of plot. Default is None
    show_legend : Optional[bool]
       Legend (color <-> cell cycle fraction) will be displayed if True. Default is False

    Returns
    -------

    """
    if plot_num is not None:
        grid_size = (5, 2)
        plot_num = plot_num % 10
        rel_pos = np.unravel_index(plot_num, grid_size)
        ax = plt.subplot2grid(grid_size, rel_pos)
    else:
        fig, ax = plt.subplots()
        
    data = dfi.copy()
    fr_col = [s for s in data_cols if 'fraction_dead' in s]
    if len(fr_col) >= 1:
        fr_col = fr_col[0]  
        data[fr_col] = -data[fr_col]
    g = data[data_cols].plot(kind='bar', stacked=True, ax=ax,
                             legend=show_legend)
    doses = ["%.2g" % d for d in data.concentration.unique()]
    ax.set_xticklabels(doses)
    if show_legend:
        plt.legend(bbox_to_anchor=(1, 0.7))
    ax.tick_params(labelsize=8)
    ax.tick_params(axis='x', labelrotation=0)
    ax.axhline(0, color='black', linewidth=1)
    ax.set_ylabel('fractions', fontsize=8)
    if title is not None:
        plt.title(title)


def plot(dfi, data_cols=['G1', 'S', 'G2', 'S_dropout',
                         'subG1', 'beyondG2',
                         'M', 'increase_fraction_dead'],
         figname='summary_distribution.pdf'):
    """Output multi-page pdf of cell cycle fractions of all 
    cell_line-drug pair in the input dataframe.

    Parameters
    ----------
    data : pandas.DataFrame 
          report dataframe from the 'run_cell_cycle_gating' function
    data_cols : Optional[list of str]
       Columns to be used to make stacked barplots
       Default is ['G1', 'S', 'G2', 'S_dropout', 'other', 'M', 'fraction_dead'].

    Returns
    -------
     
    """
    data = dfi.copy()
    data = data.groupby(['cell_line', 'agent', 'concentration'],
                        as_index=False)[data_cols].mean()
    cell_lines = data.cell_line.unique()
    agents = data.agent.unique()
    ca = list(itertools.product(cell_lines, agents))
    nb_plots = len(ca)
    pdf_pages = PdfPages(figname)
    nb_plots_per_page = 10
     
    for i, condition in enumerate(ca):
        if i % nb_plots_per_page == 0:
           fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
        cell_line = condition[0]
        agent = condition[1]
        title = "%s %s" % (cell_line, agent)
        ds = data[(data.agent == agent) & (data.cell_line == cell_line)].copy()
        if (i + 1) % nb_plots_per_page == 0 or (i + 1) == nb_plots:
            legend=True
        else:
            legend=False
        plot_stacked_bar(ds, data_cols=data_cols, plot_num=i, show_legend=legend, title=title)

        if (i + 1) % nb_plots_per_page == 0 or (i + 1) == nb_plots:
           plt.tight_layout()
           pdf_pages.savefig(fig)
           print("Completed analysis for %d out of %d conditions" %
                 (i+1, len(ca)))
           plt.close('all')
    pdf_pages.close()       
         
