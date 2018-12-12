import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import itertools
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def plot_stacked_bar(dfi, plot_num,
                     #ax=None,
                     data_cols=['G1', 'S', 'G2', 'S_dropout',
                                'other', 'M', 'fraction_dead'],
                     #condition_cols=['agent', 'cell_line', 'well'],
                     #figname='Stack_bar_plot.png',
                     #stacked=True,
                     title=None,
                     show_legend=False,
                     **kwargs):
    """Make stacked barplot. Figure is saved to the current working folder.
    
    Parameters
    ----------
    data : pd.DataFrame 
          report dataframe from the 'run_cell_cycle_gating' function
    ax : matplotlib.pyplot.Axes object 
       the axes to plot onto.
    data_cols, list of str, columns to be plotted. By default it is ['G1', 'S', 'G2', 'S_dropout', 'other', 'M'].
    condition_cols: list or str, columns names to group the data, by default it is ['agent', 'cell_line', 'well'].
    figname: str, output figure name, if None, no figure will be saved.
    stacked: bool, weather to plot stacked bars.
    """
    #ax = ax or plt.gca()
    #if ax is None:
    #    fig, ax = plt.subplots()
    grid_size = (5, 2)
    plot_num = plot_num % 10
    rel_pos = np.unravel_index(plot_num, grid_size)
    ax = plt.subplot2grid(grid_size, rel_pos)    
        
    data = dfi.copy()
    #data = data.groupby(['cell_line', 'agent', 'concentration'],
    #                    as_index=False)[data_cols].mean()
    data['fraction_dead'] = -data['fraction_dead']
    g = data[data_cols].plot(kind='bar', stacked=True, ax=ax, **kwargs,
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
                          'other', 'M', 'fraction_dead']):
     data = dfi.copy()
     data = data.groupby(['cell_line', 'agent', 'concentration'],
                        as_index=False)[data_cols].mean()
     cell_lines = data.cell_line.unique()
     agents = data.agent.unique()
     ca = list(itertools.product(cell_lines, agents))
     nb_plots = len(ca)
     pdf_pages = PdfPages('summary_distributions.pdf')
     nb_plots_per_page = 10
     
     for i, condition in enumerate(ca):
         if i % nb_plots_per_page == 0:
            fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
         cell_line = condition[0]
         agent = condition[1]
         title = "%s %s" % (cell_line, agent)
         ds = data[(data.agent == agent) & (data.cell_line == cell_line)].copy()
         if i+1 == nb_plots:
             legend=True
         else:
             legend=False
         plot_stacked_bar(ds, plot_num=i, show_legend=legend, title=title)

         if (i + 1) % nb_plots_per_page == 0 or (i + 1) == nb_plots:
            plt.tight_layout()
            pdf_pages.savefig(fig)
            print("Completed analysis for %d out of %d conditions" %
                  (i+1, len(ca)))
            plt.close('all')
     pdf_pages.close()       
         
