import pandas as pd
import xarray as xr
import numpy as np
from scipy import stats
import time
import itertools
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import datarail.utils.plate_fcts as pltfct
import datarail.utils.drug_treatment as drgtrt
from sklearn import mixture
import datarail.experimental_design.edge_fingerprint as fgprt


def Plate_bias(xray, variable='cell_count', filename=None):
    xrcc = xray[variable]

    if filename is not None:
        pdf = matplotlib.backends.backend_pdf.PdfPages(filename)

    for b in xrcc['barcode'].values:
        fig = 99
        plt.figure(fig, figsize=[10,5])
        plt.clf()

        # whole plate
        ax = plt.axes([.05, .15, .56, .7])
        himg = xrcc.loc[b].plot.imshow()
        pltfct.axis_plate(himg, xray.plate_dims)
        ax.grid(False)

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
        h.axes.axes.set_xlim([.5, xray.plate_dims[1]+.5])
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
        h.axes.axes.set_xlim([.5, xray.plate_dims[0]+.5])
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
                          xrcc.loc[b][:, np.any([xrcc.loc[b].column==i+1, # shift of 1 to avoid counting corners twice
                                                 xrcc.loc[b].column==xray.plate_dims[1]-i+1-1],axis=0)].values)
            if len(v)>0:
                m = np.append(m, v.mean())
                s = np.append(s, v.std())
                p = np.append(p, stats.ttest_1samp(v, plate_mean).pvalue)
            else:
                m = np.append(p, np.nan)
                s = np.append(p, np.nan)
                p = np.append(p, 1)


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

    vars = list(set(df.columns) - drgtrt.default_confounder_variables)
    dctrl = df.loc[df.role == 'negative_control', ['barcode'] + vars]

    # heuristic approach to define keys
    keyvars = list(set(vars) - drgtrt.possible_model_variables - \
                   {'role', 'agent', 'concentration'})
    keyvals = dctrl.loc[:,keyvars].drop_duplicates()

    plt.figure(fig, figsize=[1+2*len(keyvals),3.5])
    plt.clf()

    for ik in range(len(keyvals)):
        dkey = dctrl.loc[(dctrl.loc[:,keyvars] == keyvals.iloc[ik,:]).values.all(axis=1),:]
        dgrp = dkey.groupby('barcode')[variable]
        barcodes = dgrp.groups.keys()
        p = np.ones(len(barcodes))
        for i,b in enumerate(barcodes):
            p[i] = stats.ttest_ind(np.array(dgrp.get_group(b).values),
                                   np.concatenate([dgrp.get_group(k).values for k in dgrp.groups.keys() if k != b])).pvalue

        h = plt.axes([.1+.04/len(keyvals)+ik*.9/len(keyvals), .2, .56/len(keyvals), .7])

        m = dgrp.mean()
        s = dgrp.std()
        plt.bar(np.array(range(len(dgrp)))-.4, m)
        plt.errorbar(range(len(dgrp)), m, s, fmt='.k')
        for  i in (p<.05).nonzero()[0]:
            plt.text(i, m[i]+1.2*s[i], '*'*(1+(p[i]<.01)), horizontalalignment='center')

        h.set_title(' - '.join(keyvals.iloc[ik,:].values.astype('str')))
        h.axes.axes.set_xticks(np.array(range(len(dgrp)))-.3)
        h.axes.axes.set_xticklabels(barcodes, rotation=45)
        if ik==0: h.axes.axes.set_ylabel(variable + ' negative_control')
        h.axes.axes.set_xlim([-.5, len(dgrp)-.5])

    if filename is not None:
        pdf.savefig( fig )
        pdf.close()




def read_fingerprint(xray):
    xrcc = xray['cell_count']

    fingerprint = []
    all_pos_wells = []
    for b in xrcc['barcode'].values:
        v = np.append(xrcc.loc[b][np.any([xrcc.loc[b].row==1,
                                          xrcc.loc[b].row==xray.plate_dims[0]],axis=0),:].values,
                      xrcc.loc[b][:, np.any([xrcc.loc[b].column==2,
                                             xrcc.loc[b].column==xray.plate_dims[1]-1],axis=0)].values)

        g = mixture.GMM(n_components=2)
        g.fit([[i] for i in v])
        idx_pos = (g.means_[0]>g.means_[1])

        # get all the edge wells and their probability of being positive
        fingerprint_wells = []
        for i in range(1,xray.plate_dims[0]+1):
            fingerprint_wells.append([
                xray['well'].loc[b,i,1].values.tostring(),
                g.predict_proba([[xrcc.loc[b,i,1].values]])[0][np.argmin(g.means_)]])
            fingerprint_wells.append([
                xray['well'].loc[b,i,xray.plate_dims[1]].values.tostring(),
                g.predict_proba([[xrcc.loc[b,i,xray.plate_dims[1]].values]])[0][np.argmin(g.means_)]])

        for j in range(2,xray.plate_dims[1]):
            fingerprint_wells.append([
                xray['well'].loc[b,1,j].values.tostring(),
                g.predict_proba([[xrcc.loc[b,1,j].values]])[0][np.argmin(g.means_)]])
            fingerprint_wells.append([
                xray['well'].loc[b,xray.plate_dims[0],j].values.tostring(),
                g.predict_proba([[xrcc.loc[b,xray.plate_dims[0],j].values]])[0][np.argmin(g.means_)]])

        # decode based on the identity of the positive wells
        pos_wells = [w[0] for w in fingerprint_wells if w[1]>.5]
        fingerprint.append([b, fgprt.decode_fingerprint(pos_wells, xray.plate_dims)])
        all_pos_wells.append([b, pos_wells])

    for b in fingerprint:
        if not b[1][1]:
            print('Fingerprint failed for plate %s: read as %s !' % (b[0], b[1][0]))

    return fingerprint, all_pos_wells
