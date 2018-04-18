import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages


def plot_dose_response(df_grvalues, df_grmetrics,
                       figname='dose_response.pdf'):
    plt.ioff()
    timepoints = df_grvalues.timepoint.unique()
    pdf_pages = PdfPages(figname)
    for tp in timepoints:
        dfgrv = df_grvalues[df_grvalues.timepoint == tp]
        dfgrm = df_grmetrics[df_grmetrics.timepoint == tp]
        g = plot_dr(dfgrv, dfgrm)
        pdf_pages.savefig(g.fig)
    pdf_pages.close()


def plot_dr(df_grvalues, df_grmetric):
    agents = df_grvalues.agent.unique()
    col_wrap = np.min((5,
                       np.max((3, int(round(np.sqrt(len(agents))))))
                       ))
    num_rows = np.ceil(len(agents)/col_wrap)
    subplot_height = np.min((2, round(16 / num_rows)))
    hue_order = df_grvalues.cell_line.unique()[::-1]

    g = sns.FacetGrid(df_grvalues, col='agent', hue='cell_line',
                      # row='timepoint',
                      hue_order=hue_order, palette='husl',
                      col_wrap=col_wrap,
                      size=subplot_height, aspect=1.5,  # margin_titles=True,
                      sharey=True, sharex=False)
    timepoint = df_grvalues.timepoint.unique()[0]
    g.set(xscale='log')
    g = (g.map(plt.plot, "concentration", "GRvalue", marker='.'))
    g.set_titles(col_template='{col_name}')
    g.set_axis_labels(x_var='Concentration (uM)', y_var='GR value')
    # g.add_legend()
    labels = hue_order
    colors = sns.color_palette("husl", len(labels)).as_hex()
    handles = [patches.Patch(color=col, label=lab) for col, lab
               in zip(colors, labels)]
    xticks = [1e-3, 1e-1, 1e+1]
    for ax in g.axes.reshape(-1):
        xlim = ax.get_xlim()
        ylim_max = ax.get_ylim()[1]
        ax.plot(xlim, [0, 0])
        ax.set_xticks(xticks)
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.set_ylim((-0.6, ylim_max))
        plot_gr50(df_grmetric, ax, labels, colors)
        # plot_grmax(df_grmetric, ax, labels, colors)
    plt.subplots_adjust(hspace=0.7, wspace=0.3, bottom=0.1,
                        left=0.1, right=0.6)
    plt.legend(handles=handles, title='cell line (%s hrs)' % timepoint,
               loc='center right',
               bbox_to_anchor=(col_wrap, 1.2))
    return g


def plot_gr50(dfgr, ax, labels, colors):
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
