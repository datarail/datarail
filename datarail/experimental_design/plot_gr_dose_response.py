import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import pandas as pd


def plot_dose_response(df_grvalues, df_grmetrics=None, gr_value='GRvalue',
                       time_col='timepoint', errbar=None,
                       figname='dose_response.pdf'):
    """Plots dose response summary figures, one timepoint per pdf page.
    
    Parameters
    ----------
    df_grvalues : pandas dataframe
       GR values for each condition
    df_grmetrics : Optional[pandas dataframe]
       GR metrics summary for each condition.
    figname : Optional[str]
       Name of pdf file for saving the output.
    errbar : Optional[str]
        Default is None. `sd` to plot stdev errorbar 
    """
    plt.ioff()
    timepoints = df_grvalues.timepoint.unique()
    pdf_pages = PdfPages(figname)
    for tp in timepoints:
        dfgrv = df_grvalues[df_grvalues.timepoint == tp]
        if df_grmetrics is not None:
            dfgrm = df_grmetrics[df_grmetrics.timepoint == tp]
        else:
            dfgrm = None
        g = plot_dr(dfgrv, dfgrm, time_col, gr_value, errbar)
        pdf_pages.savefig(g.fig)
    pdf_pages.close()


def plot_dr(data, df_grmetric=None, time_col='timepoint', gr_value='GRvalue',
            errbar=None, size_override=None):
    """Plots dose response summary figure for each timepoint.

    Parameters
    ----------
    df_grvalues : pandas dataframe
       GR values for each condition
    df_grmetrics : Optional[pandas dataframe]
       GR metrics summary for each condition
    time_col : Optional[str]
      name of column containing information on duration of drug treatment
    gr_value : Optional[str]
      name of GR column to plot. Default is GRvalue. 
      Alternates include GR_static and GR_toxic.
    errbar : Optional[str]
      style in which errorbar is displayed. Default is bar
    size_override: Optional
       override default aspect ratio specified based on number of subplots.

    Returns
    -------
    g : seaborn figure object

    """
    df_grvalues = data.copy()
    df_grvalues[gr_value] = df_grvalues[gr_value].astype(float)
    dfg1 = df_grvalues.groupby(['cell_line', 'agent', 'concentration', 'timepoint'],
                                 as_index=True)[gr_value].mean().copy()
    dfg2 = df_grvalues.groupby(['cell_line', 'agent', 'concentration', 'timepoint'],
                                 as_index=True)[gr_value].std().copy()
    df_grvalues = pd.concat([dfg1, dfg2], axis=1)
    df_grvalues.columns = [gr_value, 'sd']
    df_grvalues = df_grvalues.reset_index()
    agents = df_grvalues.agent.unique()
    col_wrap = np.min((5,
                       np.max((3, int(round(np.sqrt(len(agents))))))
                       ))
    num_rows = np.ceil(len(agents)/col_wrap)
    if size_override is None:
        subplot_height = np.min((2, round(16 / num_rows)))
    else:
        subplot_height = size_override
    hue_order = df_grvalues.cell_line.unique()[::-1]

    g = sns.FacetGrid(df_grvalues, col='agent', hue='cell_line',
                      # row='timepoint',
                      hue_order=hue_order, palette='husl',
                      col_wrap=col_wrap,
                      height=subplot_height, aspect=1.5,  # margin_titles=True,
                      sharey=True, sharex=False)
    timepoint = df_grvalues[time_col].unique()[0]
    g.set(xscale='log')
    if errbar is not None:
        g = g.map(plt.errorbar, "concentration", gr_value, errbar, marker='.')
    else:
        g = g.map(plt.plot, "concentration", gr_value, marker='.')
    g.set_titles(col_template='{col_name}')
    ylabel = gr_value.replace('_', ' ')
    g.set_axis_labels(x_var=u'concentration (\u03bcM)', y_var=ylabel)
    # g.add_legend()
    labels = hue_order
    colors = sns.color_palette("husl", len(labels)).as_hex()
    handles = [patches.Patch(color=col, label=lab) for col, lab
               in zip(colors, labels)]
    xticks = [1e-3, 1e-1, 1e+1]
    for ax in g.axes.reshape(-1):
        xlim = ax.get_xlim()
        ylim_max = np.min((2, ax.get_ylim()[1]))
        ax.plot(xlim, [0, 0])
        ax.set_xticks(xticks)
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.set_ylim((-0.6, ylim_max))
        if df_grmetric is not None:
            plot_gr50(df_grmetric, ax, labels, colors)
        # plot_grmax(df_grmetric, ax, labels, colors)
    plt.subplots_adjust(hspace=0.7, wspace=0.3, bottom=0.1,
                        left=0.1, right=0.6)
    plt.legend(handles=handles, title='cell line (%s hrs)' % timepoint,
               loc='center right',
               bbox_to_anchor=(col_wrap, 1.2))
    return g


def plot_gr50(dfgr, ax, labels, colors):
    """Plots GR50 for each drug as vertical line on X-axis.

    Parameters
    ----------
    dfgr : pandas dataframe
       GR metrics summary for each condition
    ax : subplot figure object
    labels : list of str
        List of unique cell lines
    colors : array of str
        List of unique colors for each cell line. 
        Should be the same length as labels
    """
    drug = ax.get_title()
    # print(drug)
    ylim = ax.get_ylim()
    xlim = ax.get_xlim()
    dfgr_drug = dfgr[dfgr.agent == drug].copy()
    for line, color in zip(labels, colors):
        gr50 = dfgr_drug[dfgr_drug.cell_line == line]['GR50'].values[0]
        ax.plot([gr50, gr50], [ylim[0], (ylim[0] + 0.2)], color)
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)


def plot_grmax(dfgr, ax, labels, colors):
    """Plots GRmax for each drug as vertical line on Y-axis.

    Parameters
    ----------
    dfgr : pandas dataframe
       GR metrics summary for each condition
    ax : subplot figure object
    labels : list of str
        List of unique cell lines
    colors : annary of str
        List of unique colors for each cell line. 
        Should be the same length as labels
    """
    drug = ax.get_title()
    ylim = ax.get_ylim()
    xlim = ax.get_xlim()
    print(xlim)
    dfgr_drug = dfgr[dfgr.agent == drug].copy()
    for line, color in zip(labels, colors):
        grmax = dfgr_drug[dfgr_drug.cell_line == line]['GRmax'].values[0]
        ax.plot([25, 50], [grmax, grmax], color)
        ax.set_ylim(ylim)
        xlim = ax.get_xlim()
        ax.set_xlim(xlim)




def plot_fraction_dead(df_fd, y_col='fraction_dead', figname=None, errbar=None):
    """Plots fraction dead figures, one timepoint per pdf page.
    
    Parameters
    ----------
    df_fd : pandas dataframe
       fraction dead  values for each condition
    y_col : Optional[str]
       name of column to plot. Default is 'fraction_dead'.
       Alternatively, 'increase_fraction_dead' can be plotted.
    figname : Optional[str]
       Name of pdf file for saving the output.
    """
    if figname is None:
        figname = "%s.pdf" % y_col
    plt.ioff()
    timepoints = df_fd.timepoint.unique()
    pdf_pages = PdfPages(figname)
    for tp in timepoints:
        dft = df_fd[df_fd.timepoint == tp]
        g = plot_fd(dft, y_col, errbar)
        pdf_pages.savefig(g.fig)
    pdf_pages.close()



def plot_fd(data, y_col='fraction_dead', errbar=None):
    """Plots fraction dead summary figure for each timepoint.

    Parameters
    ----------
    df_fd : pandas dataframe
       Fraction dead values for each condition
    y_col : Optional[str]
       name of column to plot. Default is 'fraction_dead'.
       Alternatively, 'increase_fraction_dead' can be plotted
    errbar : Optional[str]
        Default is None. `sd` to plot stdev errorbar 

    Returns
    -------
    g : seaborn figure object
    """
    df_fd = data.copy()
    df_fd1 = df_fd.groupby(['cell_line', 'agent', 'concentration', 'timepoint'],
                           as_index=True)[y_col].mean().copy()
    df_fd2 = df_fd.groupby(['cell_line', 'agent', 'concentration', 'timepoint'],
                           as_index=True)[y_col].std().copy()
    df_fd = pd.concat([df_fd1, df_fd2], axis=1)
    df_fd.columns = [y_col, 'sd']
    df_fd = df_fd.reset_index()
    agents = df_fd.agent.unique()
    col_wrap = np.min((5,
                       np.max((3, int(round(np.sqrt(len(agents))))))
                       ))
    num_rows = np.ceil(len(agents)/col_wrap)
    subplot_height = np.min((2, round(16 / num_rows)))
    hue_order = df_fd.cell_line.unique()[::-1]

    g = sns.FacetGrid(df_fd, col='agent', hue='cell_line',
                      # row='timepoint',
                      hue_order=hue_order, palette='husl',
                      col_wrap=col_wrap,
                      height=subplot_height, aspect=1.5,  # margin_titles=True,
                      sharey=True, sharex=False)
    timepoint = df_fd.timepoint.unique()[0]
    g.set(xscale='log')
    if errbar is not None:
        g = (g.map(plt.errorbar, "concentration", y_col, errbar, marker='.'))
    else:
        g = (g.map(plt.plot, "concentration", y_col, marker='.'))
    g.set_titles(col_template='{col_name}')
    ylabel = y_col.replace('_', ' ')
    g.set_axis_labels(x_var=u'concentration (\u03bcM)', y_var=ylabel)
    # g.add_legend()
    labels = hue_order
    colors = sns.color_palette("husl", len(labels)).as_hex()
    handles = [patches.Patch(color=col, label=lab) for col, lab
               in zip(colors, labels)]
    xticks = [1e-3, 1e-1, 1e+1]
    for ax in g.axes.reshape(-1):
        xlim = ax.get_xlim()
        ylim_max = ax.get_ylim()[1]
        ax.set_xticks(xticks)
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.set_ylim((0, ylim_max))
    plt.subplots_adjust(hspace=0.7, wspace=0.3, bottom=0.1,
                        left=0.1, right=0.6)
    plt.legend(handles=handles, title='cell line (%s hrs)' % timepoint,
               loc='center right',
               bbox_to_anchor=(col_wrap, 1.2))
    return g
     
