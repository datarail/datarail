import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import os


def plot_drugs(Design, drugs):
    xarray_struct = Design.copy(deep=True)
    newpath = r'./OUTPUT/%s' % xarray_struct.barcode
    if not os.path.exists(newpath):
        os.makedirs(newpath)

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
        plt.title('%s: %s' % (xarray_struct.barcode, drug))
        plt.savefig('%s/%s_layout.png' % (newpath, drug))
        plt.clf()


def plot_all_drugs(Design, drugs=[]):

    plate_dims = Design.plate_dims
    if len(drugs)==0:
        drugs = [d for d in Design.keys() if 'DrugName' in Design[d].attrs.keys()]

    newpath = r'./OUTPUT/%s' % Design.barcode
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    colors = plt.cm.rainbow(np.linspace(0,1,len(drugs)))
    plate_RGB = np.zeros(plate_dims + [3])

    for i,drug in enumerate(drugs):
        nparray = Design[drug].values
        pos = np.nonzero(nparray)
        print pos
        concs = np.log10(nparray[pos])

        for j,c in enumerate(concs):
            plate_RGB[pos[0][j], pos[1][j], :] = \
                        np.minimum(np.maximum(colors[i][:3] + concs[j]*.1,0),1)

    plate_RGB[np.nonzero(~Design['treated_wells'].values)] = 1
    plt.imshow(plate_RGB, interpolation='None')
    plt.yticks(range(plate_dims[0]),
               list(Design[drug].coords['rows'].values))
    plt.xticks(np.arange(1, plate_dims[1]+1, 2))
    plt.gca().invert_yaxis()
    plt.title('%s - drugs' % Design.barcode)
    plt.savefig('%s/drugs_layout.png' % newpath)
    plt.clf()


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


def get_treatment_panel(design, treatment):
    drug_array = design.where(design.treatments ==
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
    treatments = np.unique(design.treatments.values)
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


def plot_layout(design):
    treatments = get_treatments(design)
    colors = plt.cm.rainbow(np.linspace(0, 1, len(treatments)))
    plate_dims = [16, 24]
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
    # plt.pcolor(plate_RGB)
    plt.imshow(plate_RGB, interpolation='None')
    plt.yticks(range(plate_dims[0]),
               list(design.coords['rows'].values))
    plt.xticks(range(plate_dims[1]),
               np.arange(1, plate_dims[1]+1))
    plt.savefig('drugs_layout.png')
    plt.clf()
    return plate_RGB
