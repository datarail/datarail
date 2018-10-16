import pandas as pd
import matplotlib.pyplot as plt

def add_median_std(data, custom_grouping = None, drug_col='agent', time_col='time', cell_col='cell_line', dose_col='concentration', gr_col='GRvalue'):
    '''
    Add calculate median and standard deviation from replicate GR values outputeed by the gr50 package.
    Parameters
    ========
    data: pd.DataFrame, gr50 output.
    (drug, dose, time, cell, GRvalue)_col: str, column names for metadata.
    custom_grouping: list, a list of column names that combined identifies each treatment condition.
    '''
    if custom_grouping is None:
        group_by = [drug_col,time_col,cell_col,dose_col]
    else:
        group_by = custom_grouping
    median_data = data.groupby(group_by)[gr_col].median()
    std_data = data.groupby(group_by)[gr_col].std()
    processed_data = pd.concat([median_data,std_data],axis = 1)
    processed_data.columns = ['median_GR','std_GR']
    processed_data.sort_values('median_GR',ascending=False,inplace=True)
    processed_data = processed_data.reset_index()
    processed_data['condition'] = processed_data[[drug_col, cell_col, dose_col, time_col]].apply(lambda x: '_'.join(x.astype(str)),axis = 1)
    return processed_data

def waterfall_plot(data, title, color_sep='cell_line', figname=None, figsize=(16,9), dose_col='concentration'):
    data = add_median_std(data)
    data.sort_values([dose_col,color_sep],ascending=False,inplace=True)
    fig_width,fig_height = (16,9)
    plt.figure(figsize=(fig_width,fig_height))
    width = fig_width/data.shape[0]/2
    num_colors = data[color_sep].unique().shape[0]
    for i, color in enumerate(data[color_sep].unique()):
        plot_data = data[data[color_sep]==color]
        sorted_dose_vector = sorted(data[dose_col].unique())
        c = ((i+1)/num_colors,0.5,1)
        color_table = pd.Series(sns.light_palette(c,input='hls',n_colors=len(sorted_dose_vector)),index = sorted_dose_vector)
        plt.bar(plot_data.index,plot_data.median_GR,yerr=plot_data.std_GR,width=width*2,color=color_table[plot_data[dose_col]],label=color)
    plt.title(title, size=16)
    plt.legend(fontsize=16)
    plt.ylabel('GRvalue',size=16)
    plt.xlabel('Treatment conditions (dose shown as color saturation, low dose to low saturation)',size=16)
    plt.xticks([])
    if figname is not None:
        plt.savefig(figname)
    else:
        plt.show()
    plt.close()

def waterfall_plot_panel(data, row_by='agent',col_by='time', color_sep='cell_line', dose_col='concentration', gr_col='GRvalue', plot_size=(4,3)):
    row_conds = data[row_by].unique()
    col_conds = data[col_by].unique()
    nrows = row_conds.shape[0]
    ncols = col_conds.shape[0]
    for row, row_value in enumerate(row_conds):
        for col, col_value in enumerate(col_conds):
            plot_data = data[(data[row_by]==row_value)&(data[col_by]==col_value)]
            title = '{}:{}'.format(str(row_value),str(col_value))
            waterfall_plot(plot_data, title, figsize=plot_size, dose_col=dose_col, color_sep=color_sep, figname=(title+'.png').replace(':','_'))

# def waterfall_plot_panel(data, row_by='agent',col_by='time', color_sep='cell_line', dose_col='concentration', gr_col='GRvalue', plot_size=(4,3)):
#     dose_order = sorted(data[dose_col].unique())
#     color_order = sorted(data[color_sep].unique())
#     g = sns.FacetGrid(data, row=row_by,col=col_by,hue=color_sep)
#     g.map(sns.barplot,dose_col, gr_col, 
#         order=dose_order, 
#         hue_order=color_order,
#         ci='sd',
#         dodge=True).add_legend()