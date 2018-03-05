import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches

def plot_dose_response(df_grvalues, df_grmetric):
    hue_order = df_grvalues.cell_line.unique()[::-1]
    g = sns.FacetGrid(df_grvalues, col='agent', hue='cell_line', #row='timepoint'
                      hue_order=hue_order, palette='husl',
                      col_wrap=3,
                      size=1, aspect=1,# margin_titles=True,
                      sharey =True, sharex=False)
    g.set(xscale='log')
    g = (g.map(plt.plot, "concentration", "GRvalue", marker='.'))
    g.set_titles(col_template='{col_name}', row_template='{row_name} hrs')
    g.set_axis_labels(x_var='Concentration (uM)', y_var='GR value')
    # g.add_legend()
    labels = hue_order
    colors = sns.color_palette("husl", len(labels)).as_hex()
    handles=[patches.Patch(color=col, label=lab) for col, lab
             in zip(colors, labels)]
    xticks = [1e-3, 1e-1, 1e+1]
    # return g
    for ax in g.axes.reshape(-1):
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()[0]
        ax.plot(xlim, [0,0])
        ax.set_xticks(xticks)
        ax.tick_params(axis='both', which='major', labelsize=6)
        plot_gr50(df_grmetric, ax, labels, colors)
        # plot_grmax(df_grmetric, ax, labels, colors)
    plt.subplots_adjust(hspace=0.7, wspace=0.3, bottom=0.1,
                        left=0.1, right=0.7)
    plt.legend(handles=handles, title='cell line', loc='center right',
               bbox_to_anchor=(5, 1.2))
    return plt


def plot_gr50(dfgr, ax, labels, colors):
    drug = ax.get_title()# .split(' | ')[1]
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

                          
