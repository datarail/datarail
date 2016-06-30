from get_hmslid import get_hmslid


def construct_design(drugs, cell_line_list, args):
    import xarray as xr
    import numpy as np

    plate_dims = args.plate_dims
    plate_rows = list(map(chr, range(65, 65+plate_dims[0])))
    plate_cols = range(1, plate_dims[1] + 1)
    Designs = xr.Dataset({k: (['rows', 'cols'], np.zeros(args.plate_dims))
                          for k in drugs.keys()},
                         coords={'rows': plate_rows, 'cols': plate_cols})

    Designs.attrs['Seed'] = args.Seed
    Designs.attrs['well_volume'] = args.well_volume
    Designs.attrs['plate_dims'] = plate_dims

    for drug in drugs:
        Designs[drug].attrs['DrugName'] = drug
        Designs[drug].attrs['Stock_conc'] = args.stock_conc
        Designs[drug].attrs['HMSLid'] = get_hmslid([drug])[drug]

    Designs['Perturbations'] = (('rows', 'cols'), np.zeros(plate_dims))
    Designs['Vehicle'] = (('rows', 'cols'), np.zeros(plate_dims))
    Designs['treated_wells'] = (('rows', 'cols'),
                                np.ones(plate_dims, dtype=bool))

    return Designs
