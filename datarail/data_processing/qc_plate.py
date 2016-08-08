import pandas as pd
import xarray as xr
import numpy as np
import datarail.utils.plate_fcts as pltfct
import matplotlib.pyplot as plt
import time

reload(pltfct)

df = pd.read_csv('../../tests/drug_response_data/OUTPUT/Example1_biased_results.tsv', sep='\t')

xray = pltfct.dfplate2xr(df)

def Plate_bias(xray):
    xrcc = xray['cell_count']

    plt.ion()

    for i,b in enumerate(xrcc['barcode'].values[9:]):
        plt.figure(i, figsize=[10, 5])
        plt.clf()

        # whole plate
        plt.axes([.05, .15, .56, .7])
        himg = xrcc.loc[b].plot.imshow()
        pltfct.axis_plate(himg, xray.plate_dims)

        # distribution by column
        h = plt.axes([.7, .1, .25, .2])
        s = xrcc.loc[b].std(dim='row')
        m = xrcc.loc[b].mean(dim='row')
        plt.errorbar(m.column, m.values, s.values)
        h.axes.axes.set_xlim([.5, xray.plate_dims[1]])
        h.axes.set_xticks(range(1, xray.plate_dims[1]+1,2))
        h.axes.axes.set_ylim(himg.get_clim())
        h.set_title('Column bias')

        # distribution by row
        h = plt.axes([.7, .4, .25, .2])
        s = xrcc.loc[b].std(dim='column')
        m = xrcc.loc[b].mean(dim='column')
        plt.errorbar(m.row, m.values, s.values)
        h.axes.axes.set_xlim([.5, xray.plate_dims[0]])
        h.axes.axes.set_ylim(himg.get_clim())
        h.axes.set_xticks(range(1, xray.plate_dims[0]+1,2))
        h.axes.set_xticklabels([chr(i) for i in
                                range(ord('A'), ord('A')+xray.plate_dims[0], 2)])
        h.set_title('Row bias')

        # distribution by distance from the edge
        m = np.array([])
        s = np.array([])
        for i in range(1,1+xray.plate_dims[0]/2):
            v = np.append(xrcc.loc[b][np.any([xrcc.loc[b].row==i,
                                              xrcc.loc[b].row==xray.plate_dims[0]-i+1],axis=0),:].values,
                          xrcc.loc[b][:, np.any([xrcc.loc[b].column==i,
                                                 xrcc.loc[b].column==xray.plate_dims[1]-i+1],axis=0)].values)
            m = np.append(m, v.mean())
            s = np.append(s, v.std())


        h = plt.axes([.7, .7, .25, .2])
        plt.errorbar(range(1,1+xray.plate_dims[0]/2), m, s)
        h.axes.axes.set_xlim([.5, xray.plate_dims[0]/2])
        h.axes.axes.set_ylim(himg.get_clim())
        h.set_title('Edge bias')

        ################################
        # save the image

    ##############################################
    # compile a single pdf with all the images


####################################
########################
### need to check for bias across the negative controls

###################################
