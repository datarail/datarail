import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def add_median_std(data, y_gr_col='GRmax', custom_grouping=None, drug_col='agent',
                   time_col='time', cell_col='cell_line'):
    '''
    Add calculate median and standard deviation from replicate GR values outputeed by the gr50 package.
    Parameters
    ========
    data: pd.DataFrame, gr50 output.
    (drug, dose, time, cell, GRvalue)_col: str, column names for metadata.
    custom_grouping: list, a list of column names that combined identifies each treatment condition.
    y_gr_col: str, column name for GR metric.
    (drug, time, and cell)_col: str, columns names for drug_name, treatment duration or cellline.

    Returns
    ========
    processed_data: pd.DataFrame, re-formmated data sorted by median 'y_gr_col' descendingly, indecies as ranks.
        Includes two columns 'median_GR' for median GR metric and 'std_GR' for standard deviation
        calculated from replicates.

    '''
    if custom_grouping is None:
        group_by = [drug_col, time_col, cell_col]
    else:
        group_by = custom_grouping
    median_data = data.groupby(group_by)[y_gr_col].median()
    std_data = data.groupby(group_by)[y_gr_col].std()
    processed_data = pd.concat([median_data, std_data], axis=1)
    processed_data.columns = ['median_GR', 'std_GR']
    processed_data.sort_values('median_GR', ascending=False, inplace=True)
    processed_data = processed_data.reset_index()
    processed_data.std_GR.fillna(0, inplace=True)
    return processed_data


def waterfall_plot(data, x_color_sep='cell_line', y_gr_col='GRmax', title_col='agent',
                   saturation_sep=None, figname=None, figsize=(16, 9), drug_col='agent',
                   time_col='time', cell_col='cell_line', ax=None, **kwargs):
    '''
    Make waterfall plot from gr50 output data or processed data from gr50 output.
    This plot can visulize two different types of metadata such as drug vs time,
    or drug vs cellline.Groups data based on 'x_color_sep' value, and saturation of
    each bar can be modified to reflect a secondary metadata.

    Parameters
    ========
    data: pd.DataFrame, processed data from gr50 output.
    x_color_sep: str, column name for X-axis which separate the primary condition.
    y_gr_col: str, column name for the GR metric values.
    saturation_sep: str or None, column name in the data which defined bar saturation level.
    figname: str or None, output figure name, if None, show figure instead.
    figsize: tuple of integers, width and height of figure.
    ax: matplotlib.Axes object, axes which the figure to be constructed on.
    title_col: str, column name where figure title will be assumed.
    (drug, dose, time, cell, GRvalue)_col: str, column names for metadata.
    y_gr_col: str, column name for GR metric.
    (drug, time, and cell)_col: str, columns names for drug_name, treatment duration or cellline.

    Returns
    ========

    '''
    # check if data is already padded median and std values, if not get those
    # values
    ax = ax or plt.gca()
    # setting xtick labels
    xtick_labels = data.iloc[:, :-2].apply(
        lambda x: '_'.join(x.astype(str).values), axis=1).values
    # sorting to ensure consistant label order in legend.
    sorting_order = pd.Series(
        {saturation_sep: False, drug_col: True, cell_col: False})
    sorting_cols = [x for x in [saturation_sep,
                                drug_col, cell_col] if x in data.columns]
    data = data.sort_values(
        sorting_cols, ascending=sorting_order[sorting_cols])
    fig_title = data[title_col].unique()
    fig_width, fig_height = figsize
    # axes wise size adjustment
    ax.figure.set_size_inches(fig_width, fig_height)

    # Setting handlers when no separation are set for color saturation or
    # color.
    if x_color_sep is None:
        x_color_sep = 'placeholder'
        data[x_color_sep] = ''
    if saturation_sep is None:
        saturation_sep = 'placeholder'
        data[saturation_sep] = ''
    # Plotting
    num_colors = data[x_color_sep].unique().shape[0]
    for i, color in enumerate(data[x_color_sep].unique()):
        plot_data = data[data[x_color_sep] == color]
        # get different saturation level within each color.
        saturation_vector = sorted(data[saturation_sep].unique())
        # If there are more than one saturation level, set level of separation
        # using sns.light_pallette.
        if len(saturation_vector) > 1:
            c = ((i + 1) / num_colors, 0.5, 1)
            # Setting color based on number of x_color_sep values, using HLS
            # value system.
            saturation_table = pd.Series(sns.light_palette(
                c, input='hls', n_colors=len(saturation_vector)), index=saturation_vector)
            ax.bar(plot_data.index, plot_data.median_GR, yerr=plot_data.std_GR,
                   color=saturation_table[plot_data[saturation_sep]], label=color)
        else:
            ax.bar(plot_data.index, plot_data.median_GR,
                   yerr=plot_data.std_GR, label=color)
    ax.axhline(0, color='black', linewidth=4)
    ax.set_title(fig_title[0], size=16)
    ax.legend(loc=1, fontsize=16)
    ax.set_ylabel(y_gr_col, size=16)
    # ax.set_xlabel(
    #     'Treatment conditions ({} shown as color saturation, low duration = low saturation)'.format(saturation_sep), size=16)
    ax.set_xticks(range(len(data)))
    ax.set_xticklabels(xtick_labels, size=14)
    if figname is not None:
        plt.savefig(figname)
        plt.close()
    return ax


def waterfall_plot_panel(data, x_color_sep='cell_line', y_gr_col='GRmax', row_by='agent', col_by=None,
                         saturation_sep=None, figname=None, figsize=(16, 9), drug_col='agent', time_col='time',
                         cell_col='cell_line', **kwargs):
    '''
    Make a panel of waterfall plot from gr50 output data or processed data from gr50 output.
    Panels can be separated by at most two types of metadata into rows and columns.

    Parameters
    ========
    data: pd.DataFrame, processed data from gr50 output.
    x_color_sep: str, column name for X-axis which separate the primary condition.
    y_gr_col: str, column name for the GR metric values.
    row_by: str, column name that separate waterfall plots by rows.
    col_by=None: str, column name that separate waterfall plots by columns. If None, plots will not be
        separated by columns, and bars will be stratified by the 'saturation_sep' parameter.
    saturation_sep: str or None, column name in the data which defined bar saturation level.
    figname: str or None, output figure name, if None, show figure instead.
    figsize: tuple of integers, width and height of figure.
    ax: matplotlib.Axes object, axes which the figure to be constructed on.
    title_col: str, column name where figure title will be assumed.
    (drug, dose, time, cell, GRvalue)_col: str, column names for metadata.
    custom_grouping: list, a list of column names that combined identifies each treatment condition.
    y_gr_col: str, column name for GR metric.
    (drug, time, and cell)_col: str, columns names for drug_name, treatment duration or cellline.

    Returns
    ========
    fig, matplotlib.pyplot.figure object. 
    '''

    if np.sum(['median' in x for x in data.columns]) == 0:
        data = add_median_std(data, **kwargs)
    # setting rows and cols in subplots based on row_by and col_by parameters.
    row_conds = data[row_by].unique()
    nrows = row_conds.shape[0]
    ncols = 1
    if col_by is not None:
        col_conds = data[col_by].unique()
        ncols = col_conds.shape[0]
    fig, axes = plt.subplots(nrows, ncols, sharex=False, sharey=True)
    # iterate through each axes and plot into it.
    for row, row_cond in enumerate(row_conds):
        # If additional col level separation are present.
        # Make a nrows X ncols plot.
        if ncols > 1:
            for col, col_cond in enumerate(col_conds):
                plot_ax = axes[row, col]
                plot_data = data[(data[row_by] == row_cond) &
                                 (data[col_by] == col_cond)]
                plot_data.index = np.arange(plot_data.shape[0])
                waterfall_plot(plot_data, title_col=row_by, ax=plot_ax, saturation_sep=saturation_sep,
                               x_color_sep=x_color_sep, figsize=figsize, y_gr_col=y_gr_col, **kwargs)
                # Supress x label and legends.
                plot_ax.set_xlabel('')
                plot_ax.legend(bbox_to_anchor=(1.05, 0.7), fontsize=16)
                plot_ax.set_title('{}:{}'.format(
                    str(row_cond), str(col_cond)), size=16)
                plot_ax.legend_.remove()
                plot_ax.tick_params(labelrotation=15)
        # If no col level separation are present, make a nrows X 1 plot.
        else:
            plot_ax = axes[row]
            plot_data = data[data[row_by] == row_cond]
            plot_data.index = np.arange(plot_data.shape[0])
            waterfall_plot(plot_data, title_col=row_by, ax=plot_ax, saturation_sep=saturation_sep,
                           x_color_sep=x_color_sep, figsize=figsize, y_gr_col=y_gr_col, **kwargs)
            plot_ax.set_xlabel('')
            plot_ax.legend_.remove()
            plot_ax.tick_params(labelrotation=15)
    # fixing legend and layouts
    num_legend_cols = len(plot_ax.get_legend_handles_labels()[1])
    legend_labels = plot_ax.get_legend_handles_labels()[1]
    legend_handles = plot_ax.get_legend_handles_labels()[0]
    fig.legend(legend_handles, legend_labels,
               loc='center', bbox_to_anchor=(1.1, 0.5), fontsize='x-large')
    plt.tight_layout()
    if figname is not None:
        plt.savefig(figname)
    return fig
