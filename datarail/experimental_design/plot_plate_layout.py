import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec


def plot_feature(df, plate=None, feature='cell_line', ax=None):
    if ax is None:
        ax = plt.figure()
    if plate is None:
        plate = df.barcode.unique()[1]
    dfp = df[df.barcode == plate].copy()
    dfp[feature] = dfp[feature].fillna('untreated')
    hue_order = dfp[feature].unique()[::-1]
    labels = hue_order
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
    cg.legend(loc=(1, 0.5))
    # ax.subplots_adjust(right=0.8)
    ax.set_title(plate)
    return dfpv


def plot_concentration(df, plate=None, ax=None):
    if ax is None:
        fig, ax = plt.subplot()
    if plate is None:
        plate = df.barcode.unique()[1]
    dfp = df[df.barcode == plate].copy()
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
    sns.heatmap(dfpv, norm=LogNorm(s.min(), s.max()),
                cmap=cmap,
                linewidth=0.3, linecolor='white',
                cbar_kws={"ticks": [1e-4, 1e-3, 1e-2, 1e-1,
                                    1, 10, 1e2, 1e3, 1e4, 1e5],
                          "label": 'concentration (um)'},
                vmin=1e-6, cbar=True, ax=ax)
    ax.set_title(plate)
    return dfpv


def plot_summary(df, plate=None):
    if plate is None:
        plate = df.barcode.unique()[1]
    fig = plt.figure(figsize=(20, 4))
    plt.subplots_adjust(right=0.95, left=0.05, wspace=0.5)
    gridspec.GridSpec(1, 3)
    ax1 = plt.subplot2grid((1, 3), (0, 0), colspan=1, rowspan=1)
    ax2 = plt.subplot2grid((1, 3), (0, 1), colspan=1, rowspan=1)
    ax3 = plt.subplot2grid((1, 3), (0, 2), colspan=1, rowspan=1)
    plot_feature(df, plate=plate, feature='cell_line', ax=ax1)
    plot_feature(df, plate=plate, feature='agent', ax=ax2)
    plot_concentration(df, plate=plate, ax=ax3)
    return fig
