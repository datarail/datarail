import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def plot_feature(df, plate=None, feature='cell_line', ax=None):
    """ Plots layout of agents(drug) or cell lines for a single plate
    Parameters
    ----------
    df: pandas dataframe
       metadata table with barcode, well, agent, cell_line,
       and concentration as columns
    plate: str
      identifier of plate in barcode column
    feature: str
      feature to plot ('cell_line' or 'agent')
    ax: subplot object
      provides positional reference for plot in master pdf
    Returns
    -------
    dfpv: pandas dataframe
      metadata table pivoted such that it has plate dimensions
      value in dfpv is the feature selected.
    """
    if ax is None:
        fig, ax = plt.subplots()
    if plate is None:
        plate = df.barcode.unique()[0]
    dfp = df[df.barcode == plate].copy()
    tp = dfp.timepoint.unique()[0]
    tp = tp.replace('.0', '')
    if tp != 'time0_ctrl':
        tp = "%s hrs" % tp
    if feature != 'cell_line':
        agent_cols = [a for a in dfp.columns.tolist()
                      if 'agent' in a]
        dfp[agent_cols] = dfp[agent_cols].fillna('untreated')
        labels = np.sort(dfp[agent_cols[0]].unique())
    else:
        dfp[feature] = dfp[feature].fillna('untreated')
        labels = np.sort(dfp[feature].unique())
    colrs = sns.color_palette("husl", len(labels)).as_hex()
    colrs[np.argmax(labels == 'untreated')] = 'whitesmoke'
    if 'DMSO' in labels:
        colrs[np.argmax(labels == 'DMSO')] = 'black'
    # update labels and colors
    nlabels = np.sort(dfp[feature].unique())
    colrs = np.array(colrs)
    colrs = colrs[[l in nlabels for l in labels]]
    color_dict = {lab: col for lab, col in zip(nlabels, colrs)}
    colrs = list(colrs)
    cmap = colors.ListedColormap(colrs)

    label_map = {l: i for i, l in enumerate(nlabels)}
    dfp[feature] = dfp[feature].map(label_map)
    dfp['row'] = [s[0] for s in dfp.well.tolist()]
    dfp['column'] = [s[1:] for s in dfp.well.tolist()]
    dfpv = dfp.pivot(index='row', columns='column', values=feature)
    cg = sns.heatmap(dfpv, cmap=cmap,
                     linewidth=0.3, linecolor='white',
                     cbar=False, ax=ax)
    for label in color_dict.keys():
        cg.bar(0, 0, color=color_dict[label],
               label=label, linewidth=0)
    cg.legend(loc=(0.3, -0.6))
    ax.set_title("%s (%s)" % (plate, tp))
    ax.set_xlabel('')
    ax.set_ylabel('')
    return dfpv


def plot_concentration(df, feature='concentration',
                       plate=None, ax=None):
    """ Plots layout of agent(drug) concentrations for a single plate
    Parameters
    ----------
    df: pandas dataframe
       metadata table with barcode, well, agent, cell_line, 
       and concentration as columns
    plate: str
      identifier of plate in barcode column
    feature: str
      feature to plot (default is concentration)
    ax: subplot object
      provides positional reference for plot in master pdf
    Returns
    -------
    dfpv: pandas dataframe
      metadata table pivoted such that it has plate dimensions
      value in dfpv is the feature (concentration) selected.
    """
    if ax is None:
        fig, ax = plt.subplots()
    if plate is None:
        plate = df.barcode.unique()[0]
    dfp = df[df.barcode == plate].copy()
    tp = dfp.timepoint.unique()[0]
    tp = tp.replace('.0', '')
    if tp != 'time0_ctrl':
        tp = "%s hrs" % tp
    dfp[feature] = dfp[feature].fillna(0)
    cs = np.sort(dfp[feature].unique())
    if len(cs) > 1:
        cmin = np.log10(cs[1])
        cmax = np.log10(cs[-1])
        ticks = [10 ** c for c in np.arange(cmin, cmax+1)]
        cmd = 10 ** (cmin-1)
    else:
        ticks = [0]
        cmd = 1e-6
    dfp[feature] = dfp[feature].replace([0], cmd)
    labels = dfp[feature].unique()
    colrs = sns.color_palette("YlOrRd", len(labels)).as_hex()
    colrs[np.argmax(labels == cmd)] = 'whitesmoke'
    cmap = colors.ListedColormap(colrs)
    dfp['row'] = [s[0] for s in dfp.well.tolist()]
    dfp['column'] = [s[1:] for s in dfp.well.tolist()]
    dfpv = dfp.pivot(index='row', columns='column', values=feature)
    s = dfpv.values
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    sns.heatmap(dfpv, norm=LogNorm(s.min(), s.max()),
                cmap=cmap,
                linewidth=0.3, linecolor='white',
                cbar_kws={"ticks": ticks,
                          "label": 'concentration (um)'},
                vmin=cmd, cbar=True, ax=ax, cbar_ax=cax)
    ax.set_title("%s (%s)" % (plate, tp))
    ax.set_xlabel('')
    ax.set_ylabel('')
    return dfpv


def plot_summary(df, figname='plate_layout.pdf'):
    """ Generated multi-page pdf (one page per plate) of plate layouts.
    Each page has the layout of cell_lines (column1), agents(column2)
    and corresponding concentrations (column3)
    Parameters
    ----------
    df: pandas dataframe
       metadata table with barcode, well, agent, cell_line,
       and concentration as columns
    figname: str
      filename for saving the output pdf
    """
    plt.ioff()
    plates = df.barcode.unique()
    pdf_pages = PdfPages(figname)
    for plate in plates:
        agent_cols = [a for a in df.columns.tolist()
                      if 'agent' in a]
        conc_cols = [c for c in df.columns.tolist()
                     if 'concentration' in c]
        grid_height = len(agent_cols)
        grid_dims = (grid_height,  3)
        fig = plt.figure(figsize=(16, 6 * grid_height), dpi=100)
        bt = 0.5 - 0.1 * grid_height
        plt.subplots_adjust(left=0.03, right=0.9,
                            wspace=0.2, hspace=0.8,
                            bottom=bt, top=0.95)
        gridspec.GridSpec(*grid_dims)
        ax1 = plt.subplot2grid(grid_dims, (0, 0), colspan=1, rowspan=1)
        plot_feature(df, plate=plate, feature='cell_line', ax=ax1)
        for ai, agent in enumerate(agent_cols):
            ax = plt.subplot2grid(grid_dims, (ai, 1),
                                  colspan=1, rowspan=1)
            plot_feature(df, plate=plate, feature=agent, ax=ax)
        for ci, conc in enumerate(conc_cols):
            ax = plt.subplot2grid(grid_dims, (ci, 2),
                                  colspan=1, rowspan=1)
            plot_concentration(df, plate=plate, feature=conc, ax=ax)
        pdf_pages.savefig(fig)
    pdf_pages.close()
    return
