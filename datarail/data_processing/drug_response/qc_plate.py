import pandas as pd
import xarray as xr
import numpy as np
from scipy import stats
import time
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import datarail.utils.plate_fcts as pltfct
import datarail.utils.drug_treatment as drgtrt

reload(drgtrt)

def Plate_bias(xray, variable='cell_count', filename=None):
    xrcc = xray[variable]

    if filename is not None:
        pdf = matplotlib.backends.backend_pdf.PdfPages(filename)

    for i,b in enumerate(xrcc['barcode'].values):
        fig = 99
        plt.figure(fig, figsize=[10,5])
        plt.clf()

        # whole plate
        plt.axes([.05, .15, .56, .7])
        himg = xrcc.loc[b].plot.imshow()
        pltfct.axis_plate(himg, xray.plate_dims)

        plate_mean = xrcc.loc[b].values.mean()

        # distribution by column
        h = plt.axes([.7, .1, .25, .2])
        s = xrcc.loc[b].std(dim='row')
        m = xrcc.loc[b].mean(dim='row')
        plt.errorbar(m.column, m.values, s.values)
        plt.plot(m.column[[0, -1]], 2*[plate_mean], '-k', linewidth=2)
        for i in m.column:
            p = stats.ttest_1samp(xrcc.loc[b,:,i], plate_mean).pvalue
            plt.text(i, m.loc[i]+1.2*s.loc[i], '*'*(1*(p<.05)+(p<.01)), horizontalalignment='center')
        h.axes.axes.set_xlim([.5, xray.plate_dims[1]])
        h.axes.set_xticks(range(1, xray.plate_dims[1]+1,2))
        h.axes.axes.set_ylim(himg.get_clim())
        h.set_title('Column bias')

        # distribution by row
        h = plt.axes([.7, .4, .25, .2])
        s = xrcc.loc[b].std(dim='column')
        m = xrcc.loc[b].mean(dim='column')
        plt.errorbar(m.row, m.values, s.values)
        plt.plot(m.row[[0, -1]], 2*[plate_mean], '-k', linewidth=2)
        for i in m.row:
            p = stats.ttest_1samp(xrcc.loc[b,i,:], plate_mean).pvalue
            plt.text(i, m.loc[i]+1.2*s.loc[i], '*'*(1*(p<.05)+(p<.01)), horizontalalignment='center')
        h.axes.axes.set_xlim([.5, xray.plate_dims[0]])
        h.axes.axes.set_ylim(himg.get_clim())
        h.axes.set_xticks(range(1, xray.plate_dims[0]+1,2))
        h.axes.set_xticklabels([chr(i) for i in
                                range(ord('A'), ord('A')+xray.plate_dims[0], 2)])
        h.set_title('Row bias')

        # distribution by distance from the edge
        m,s,p = 3*[np.array([])]
        for i in range(1,1+xray.plate_dims[0]/2):
            v = np.append(xrcc.loc[b][np.any([xrcc.loc[b].row==i,
                                              xrcc.loc[b].row==xray.plate_dims[0]-i+1],axis=0),:].values,
                          xrcc.loc[b][:, np.any([xrcc.loc[b].column==i,
                                                 xrcc.loc[b].column==xray.plate_dims[1]-i+1],axis=0)].values)
            m = np.append(m, v.mean())
            s = np.append(s, v.std())
            p = np.append(p, stats.ttest_1samp(v, plate_mean).pvalue)


        h = plt.axes([.7, .7, .25, .2])
        plt.errorbar(range(1,1+xray.plate_dims[0]/2), m, s)
        plt.plot([1, xray.plate_dims[0]/2], 2*[plate_mean], '-k', linewidth=2)
        for i in (p<.05).nonzero()[0]:
            plt.text(i+1, m[i]+1.2*s[i], '*'*(1+(p[i]<.01)), horizontalalignment='center')
        h.axes.axes.set_xlim([.5, .5+xray.plate_dims[0]/2])
        h.axes.axes.set_ylim(himg.get_clim())
        h.set_title('Edge bias')

        # save the image
        if filename is not None:
            pdf.savefig( fig )

    if filename is not None:
        pdf.close()



### need to check for bias across the negative controls

def Negative_control_bias(df, variable='cell_count', filename=None):

    if filename is not None:
        pdf = matplotlib.backends.backend_pdf.PdfPages(filename)
    fig = 199

    vars = set(df.columns) - (drgtrt.default_confounder_varaibles - {'barcode'})
    dctrl = df.loc[df.role == 'negative_control', vars]

    dgrp = dctrl.groupby(['cell_line', 'barcode'])[variable]

    plt.figure(fig, figsize=[8,5])
    plt.clf()
    h = plt.axes([.12, .2, .8, .75])

    m = dgrp.mean();
    plt.bar(np.array(range(len(dgrp)))-.4, m)
    plt.errorbar(range(len(dgrp)), m, dgrp.std(), fmt='.k')

    xlabels = [m.index.levels[0][i] + ' ' + m.index.levels[1][j] for i,j in
               zip(m.index.labels[0], m.index.labels[1])]

    h.axes.axes.set_xticks(np.array(range(len(dgrp)))-.3)
    h.axes.axes.set_xticklabels(xlabels, rotation=45)
    h.axes.axes.set_ylabel(variable + ' negative_control')

    if filename is not None:
        pdf.savefig( fig )
        pdf.close()
