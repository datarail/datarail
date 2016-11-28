import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import os


def get_treatment_panel(design, treatment):
    drug_array = design.where(design.agents ==
                              treatment).concentrations.values
    plate_dims = [drug_array.shape[0], drug_array.shape[1]]
    drug_array_flat = np.zeros(plate_dims)
    combo_k = drug_array.shape[2]
    for k in range(combo_k):
        a = drug_array[:, :, k]
        a[np.isnan(a)] = 0
        drug_array_flat += a
    return drug_array_flat


def get_treatments(design):
    treatments = np.unique(design.agents.values)
    treatments_list = [tr for tr in treatments if tr != '']
    return treatments_list


def get_treated_panel(design):
    treated_array = design.concentrations.values
    plate_dims = [treated_array.shape[0], treated_array.shape[1]]
    treated_panel_flat = np.zeros(plate_dims)
    combo_k = treated_array.shape[2]
    for k in range(combo_k):
        a = treated_array[:, :, k]
        a[np.isnan(a)] = 0
        treated_panel_flat += a
    return treated_panel_flat


def plot_layout(design, filename=None):
    treatments = get_treatments(design)
    colors = plt.cm.rainbow(np.linspace(0, 1, len(treatments)))
    plate_dims = [design.dims['rows'], design.dims['columns']]
    plate_RGB = np.zeros(plate_dims + [3])

    for i, treatment in enumerate(treatments):
        treatment_array = get_treatment_panel(design, treatment)
        pos = np.nonzero(treatment_array)
        concs = np.log10(treatment_array[pos])
        for j, c in enumerate(concs):
            plate_RGB[pos[0][j], pos[1][j], :] = \
                        np.minimum(np.maximum(colors[i][:3] +
                                              concs[j]*.1, 0), 1)

    treated_panel = get_treated_panel(design)
    untreated_wells = np.nonzero(treated_panel == 0)
    plate_RGB[untreated_wells] = 1
    plt.imshow(plate_RGB, interpolation='None')
    plt.yticks(range(plate_dims[0]),
               list(design.coords['rows'].values))
    plt.xticks(range(plate_dims[1]),
               np.arange(1, plate_dims[1]+1))
    plt.title('layout of %s' % design.plates.values)
    if not filename:
        plt.show()
    elif filename:
        plt.savefig(filename)


def plot_drug(design, drug, filename=None):
    drug_array = design.where(design.agents ==
                              drug).concentrations[:, :, 0].values
    vmax = np.nanmax(drug_array)*1.2
    drug_array[np.isnan(drug_array)] = 1e6
    plate_dims = drug_array.shape
    panel = xr.DataArray(drug_array,
                         coords={'rows': range(plate_dims[0]),
                                 'columns': range(1, plate_dims[1]+1)})
    panel.plot(cmap='hot', edgecolor='k',
               norm=LogNorm(vmin=drug_array.min()/1.1, vmax=vmax),
               vmin=drug_array.min()/1.1, vmax=vmax)
    plt.yticks(range(plate_dims[0]),
               list(design.coords['rows'].values))
    plt.xticks(np.arange(1, plate_dims[1]+1, 2))
    plt.gca().invert_yaxis()
    plt.title('layout of %s (plate %s)' % (drug,design.plates.values))
    if not filename:
        plt.show()
    elif filename:
        plt.savefig(filename)



def plot_control_wells(xarray_struct):
    newpath = r'OUTPUT/%s' % xarray_struct.barcode
    if not os.path.exists(newpath):
        os.makedirs(newpath)

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
    plt.savefig('%s/control_wells_layout.png' % newpath)
    plt.clf()
