import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages


def plot_feature(df, plate=None, feature='cell_line', ax=None):
    if ax is None:
        ax = plt.figure()
    if plate is None:
        plate = df.barcode.unique()[1]
    dfp = df[df.barcode == plate].copy()
    tp = dfp.timepoint.unique()[0]
    tp = tp.replace('.0', ' hrs')
    dfp[feature] = dfp[feature].fillna('untreated')
    labels = np.sort(dfp[feature].unique())
    colrs = sns.color_palette("husl", len(labels)).as_hex()
    colrs[np.argmax(labels == 'untreated')] = 'whitesmoke'
    color_dict = {lab: col for lab, col in zip(labels, colrs)}

    cmap = colors.ListedColormap(colrs)

    label_map = {l: i for i, l in enumerate(labels)}
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
    cg.legend(loc=(0.3, -0.7))
    ax.set_title("%s (%s)" % (plate, tp))
    ax.set_xlabel('')
    ax.set_ylabel('')
    return dfpv


def plot_concentration(df, plate=None, ax=None):
    if ax is None:
        fig, ax = plt.subplots()
    if plate is None:
        plate = df.barcode.unique()[1]
    dfp = df[df.barcode == plate].copy()
    tp = dfp.timepoint.unique()[0]
    tp = tp.replace('.0', ' hrs')
    dfp['concentration'] = dfp['concentration'].fillna(0)
    dfp['concentration'] = dfp['concentration'].replace([0], 1e-6)
    labels = dfp.concentration.unique()
    colrs = sns.color_palette("YlOrRd", len(labels)).as_hex()
    colrs[np.argmax(labels == 1e-6)] = 'whitesmoke'
    cmap = colors.ListedColormap(colrs)
    dfp['row'] = [s[0] for s in dfp.well.tolist()]
    dfp['column'] = [s[1:] for s in dfp.well.tolist()]
    dfpv = dfp.pivot(index='row', columns='column', values='concentration')
    s = dfpv.values
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    sns.heatmap(dfpv, norm=LogNorm(s.min(), s.max()),
                cmap=cmap,
                linewidth=0.3, linecolor='white',
                cbar_kws={"ticks": [1e-4, 1e-3, 1e-2, 1e-1,
                                    1, 10, 1e2, 1e3, 1e4, 1e5],
                          "label": 'concentration (um)'},
                vmin=1e-6, cbar=True, ax=ax, cbar_ax=cax)
    ax.set_title("%s (%s)" % (plate, tp))
    ax.set_xlabel('')
    ax.set_ylabel('')
    return dfpv


def plot_summary(df):
    plt.ioff()
    plates = df.barcode.unique()
    pdf_pages = PdfPages('plate_layout.pdf')
    for plate in plates:
        fig = plt.figure(figsize=(16, 6), dpi=100)
        plt.subplots_adjust(left=0.03, right=0.95,
                            wspace=0.2, hspace=0.2,
                            bottom=0.4, top=0.88)
        gridspec.GridSpec(1, 3)
        ax1 = plt.subplot2grid((1, 3), (0, 0), colspan=1, rowspan=1)
        ax2 = plt.subplot2grid((1, 3), (0, 1), colspan=1, rowspan=1)
        ax3 = plt.subplot2grid((1, 3), (0, 2), colspan=1, rowspan=1)
        plot_feature(df, plate=plate, feature='cell_line', ax=ax1)
        plot_feature(df, plate=plate, feature='agent', ax=ax2)
        plot_concentration(df, plate=plate, ax=ax3)
        pdf_pages.savefig(fig)
    pdf_pages.close()
    return
