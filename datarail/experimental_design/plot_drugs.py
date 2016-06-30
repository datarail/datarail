import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def plot(xarray_struct, drugs):
    for drug in drugs:
        nparray = xarray_struct[drug].values
        nparray[nparray == 0] = 1e-2
        panel = xr.DataArray(nparray,
                             coords={'rows': range(16), 'cols': range(24)})

        panel.plot(cmap='jet', edgecolor='k',
                   norm=LogNorm(vmin=nparray.min(), vmax=nparray.max()))
        plt.yticks(range(16), list(xarray_struct[drug].coords['rows'].values))
        plt.savefig('layout_images/%s_layout.png' % drug)
        plt.clf()
