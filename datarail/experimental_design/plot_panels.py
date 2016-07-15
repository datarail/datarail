import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


def plot_drugs(xarray_struct, drugs, rep):
    for drug in drugs:
        nparray = xarray_struct[drug].values
        plate_dims = nparray.shape
        nparray[nparray == 0] = 1e-5
        panel = xr.DataArray(nparray,
                             coords={'rows': range(plate_dims[0]),
                                     'cols': range(1, plate_dims[1]+1)})

        panel.plot(cmap='jet', edgecolor='k',
                   norm=LogNorm(vmin=nparray.min(), vmax=nparray.max()))
        plt.yticks(range(plate_dims[0]),
                   list(xarray_struct[drug].coords['rows'].values))
        plt.xticks(np.arange(1, plate_dims[1]+1, 2))
        plt.gca().invert_yaxis()
        plt.title('%s_rep_%d' % (drug, rep))
        plt.savefig('ld_layout_images/%s_layout_rep%d.png' % (drug, rep))
        plt.clf()


def plot_control_wells(xarray_struct):
    bool_array = xarray_struct['control_wells'].values
    plate_dims = bool_array.shape
    panel = xr.DataArray(bool_array,
                         coords={'rows': range(plate_dims[0]),
                                 'cols': range(1, plate_dims[1]+1)})
    panel.plot(cmap='jet', edgecolor='k')
    plt.yticks(range(plate_dims[0]),
               list(xarray_struct['control_wells'].coords['rows'].values))
    plt.xticks(np.arange(1, plate_dims[1]+1, 2))
    plt.gca().invert_yaxis()
    plt.title('Control wells')
    plt.savefig('ld_layout_images/control_wells_layout.png')
    plt.clf()
