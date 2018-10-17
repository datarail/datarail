import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def add_median_std(data, custom_grouping = None, drug_col='agent', time_col='time', cell_col='cell_line', gr_col='GRmax'):
    '''
    Add calculate median and standard deviation from replicate GR values outputeed by the gr50 package.
    Parameters
    ========
    data: pd.DataFrame, gr50 output.
    (drug, dose, time, cell, GRvalue)_col: str, column names for metadata.
    custom_grouping: list, a list of column names that combined identifies each treatment condition.
    '''
    if custom_grouping is None:
        group_by = [drug_col,time_col,cell_col]
    else:
        group_by = custom_grouping
    median_data = data.groupby(group_by)[gr_col].median()
    std_data = data.groupby(group_by)[gr_col].std()
    processed_data = pd.concat([median_data,std_data],axis = 1)
    processed_data.columns = ['median_GR','std_GR']
    processed_data.sort_values('median_GR',ascending=False,inplace=True)
    processed_data = processed_data.reset_index()
    processed_data.std_GR.fillna(0,inplace=True)
    return processed_data

def waterfall_plot(data, singular_cond_col='agent', ax=None, saturation_sep='time', color_sep='cell_line', figname=None, figsize=(16,9), gr_col='GRmax', drug_col='agent', time_col='time', cell_col='cell_line', **kwargs): 
    # check if data is already padded median and std values, if not get those values
    ax = ax or plt.gca()
    if np.sum(['median' in x for x in data.columns]) == 0:
        data = add_median_std(data,**kwargs)
    # sorting to ensure consistant label order in legend.
    sorting_cols = [x for x in [saturation_sep,drug_col,cell_col] if x in data.columns]
    data = data.sort_values(sorting_cols,ascending=False)

    fig_title = data[singular_cond_col].unique()
    fig_width, fig_height = figsize
    # axes wise size adjustment
    ax.figure.set_size_inches(fig_width, fig_height)

    # Setting handlers when no separation are set for color saturation or color. 
    if color_sep is None:
        color_sep = 'placeholder'
        data[color_sep] = ''
    if saturation_sep is None:
        saturation_sep = 'placeholder'
        data[saturation_sep] = ''
    # Plotting
    num_colors = data[color_sep].unique().shape[0]
    for i, color in enumerate(data[color_sep].unique()):
        plot_data = data[data[color_sep]==color]
        # get different saturation level within each color.
        saturation_vector = sorted(data[saturation_sep].unique())
        # If there are more than one saturation level, set level of separation using sns.light_pallette.
        if len(saturation_vector)>1:
            c = ((i+1)/num_colors,0.5,1)
            # Setting color based on number of color_sep values, using HLS value system.
            saturation_table = pd.Series(sns.light_palette(c,input='hls',n_colors=len(saturation_vector)),index = saturation_vector)
            ax.bar(plot_data.index,plot_data.median_GR,yerr=plot_data.std_GR,color=saturation_table[plot_data[saturation_sep]],label=color)
        else:
            ax.bar(plot_data.index,plot_data.median_GR,yerr=plot_data.std_GR,label=color)

    ax.set_title(fig_title[0], size=16)
    ax.legend(loc=1, fontsize=16)
    ax.set_ylabel(gr_col,size=16)
    ax.set_xlabel('Treatment conditions (treatment time shown as color saturation, low duration = low saturation)',size=16)
    ax.set_xticks([])
    if figname is not None:
        plt.savefig(figname)
    # plt.close()

def waterfall_plot_panel(data, row_by='agent',col_by=None, saturation_sep='time', color_sep='cell_line', figname=None, figsize=(16,9), gr_col='GRmax', drug_col='agent', time_col='time', cell_col='cell_line', **kwargs):

    if np.sum(['median' in x for x in data.columns]) == 0:
        data = add_median_std(data,**kwargs)
    # setting rows and cols in subplots based on row_by and col_by parameters.
    row_conds = data[row_by].unique()
    nrows = row_conds.shape[0]
    ncols = 1
    if col_by is not None:
        col_conds = data[col_by].unique()
        ncols = col_conds.shape[0]
    fig,axes = plt.subplots(nrows, ncols, sharex=True, sharey=True)
    # iterate through each axes and plot into it. 
    for row, row_cond in enumerate(row_conds):
        # If additional col level separation are present. Make a nrows X ncols plot. 
        if ncols>1:
            for col, col_cond in enumerate(col_conds):
                plot_ax = axes[row,col]
                plot_data = data[(data[row_by]==row_cond)&(data[col_by]==col_cond)]
                plot_data.index = np.arange(plot_data.shape[0])
                waterfall_plot(plot_data, singular_cond_col=row_by, ax=plot_ax, saturation_sep=saturation_sep, color_sep=color_sep, figsize=figsize, **kwargs)
                # Supress x label and legends.
                plot_ax.set_xlabel('')
                plot_ax.legend(bbox_to_anchor=(1.05,0.7), fontsize=16)
                plot_ax.set_title('{}:{}'.format(str(row_cond),str(col_cond)), size=16)
        # If no col level separation are present, make a nrows X 1 plot.
        else:
            plot_ax = axes[row]
            plot_data = data[data[row_by]==row_cond]
            plot_data.index = np.arange(plot_data.shape[0])
            waterfall_plot(plot_data, singular_cond_col=row_by, ax=plot_ax, saturation_sep=saturation_sep, color_sep=color_sep, figsize=figsize, **kwargs)
            plot_ax.set_xlabel('')
    plt.tight_layout()
    if figname is not None:
        plt.savefig(figname)
    return fig