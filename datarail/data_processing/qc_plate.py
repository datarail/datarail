import pandas as pd
import xarray as xr
import datarail.utils.plate_fcts as pltfct
import matplotlib.pyplot as plt
import time

df = pd.read_csv('../../tests/drug_response_data/OUTPUT/Example1_biased_results.tsv', sep='\t')

xray = pltfct.dfplate2xr(df)

xrcc = xray['cell_count']

plt.ion()

for i,b in enumerate(xrcc['barcode'].values[11:]):
    plt.figure(i, figsize=[10, 5])
    plt.clf()

    plt.axes([.07, .15, .6, .7])
    xrcc.loc[b].plot.imshow()

    plt.axes([.7, .1, .25, .2])
    (xrcc.loc[b].std(dim='row')/xrcc.loc[b].mean(dim='row')).plot()

    plt.axes([.7, .4, .25, .2])
    (xrcc.loc[b].std(dim='column')/xrcc.loc[b].mean(dim='column')).plot()
