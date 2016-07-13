import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def plot_drugs(xarray_struct, drugs, rep):
    for drug in drugs:
        nparray = xarray_struct[drug].values
        nparray[nparray == 0] = 1e-5
        panel = xr.DataArray(nparray,
                             coords={'rows': range(16), 'cols': range(24)})

        panel.plot(cmap='jet', edgecolor='k',
                   norm=LogNorm(vmin=nparray.min(), vmax=nparray.max()))
        plt.yticks(range(16), list(xarray_struct[drug].coords['rows'].values))
        plt.savefig('ld_layout_images/%s_layout_rep%d.png' % (drug, rep))
        plt.clf()


def plot_control_wells(xarray_struct):
    bool_array = xarray_struct['control_wells'].values
    panel = xr.DataArray(bool_array,
                         coords={'rows': range(16), 'cols': range(24)})
    panel.plot(cmap='jet', edgecolor='k')
    plt.yticks(range(16),
               list(xarray_struct['control_wells'].coords['rows'].values))
    plt.savefig('ld_layout_images/control_wells_layout.png')
    plt.clf()
