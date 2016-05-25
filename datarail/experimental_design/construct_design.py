def construct_design(drug_list, cell_line_list, args):
    import xarray as xr
    import numpy as np

    plate_dims = args.plate_dims
    plate_rows = list(map(chr, range(65, 65+plate_dims[0])))
    plate_cols = range(1, plate_dims[1] + 1)

    
    Designs = xr.Dataset({k.name: (['rows', 'cols'], np.zeros(plate_dims))
                          for k in drug_list},
                         coords={'rows': plate_rows, 'cols': plate_cols})


    # Designs = xr.Dataset({k: (['rows', 'cols'], np.zeros(args.plate_dims))
    #                       for k in drugs.keys()},
    #                      coords={'rows': plate_rows, 'cols': plate_cols})

    Designs.attrs['Seed'] = args.Seed
    Designs.attrs['well_volume'] = args.well_volume
    Designs.attrs['plate_dims'] = plate_dims

    for drug in drug_list:
        Designs[drug.name].attrs['DrugName'] = drug.name
        Designs[drug.name].attrs['Stock_conc'] = drug.stock_concentration
        Designs[drug.name].attrs['HMSLid'] = drug.database_id

    Designs['Perturbations'] = (('rows', 'cols'), np.zeros(plate_dims))
    Designs['Vehicle'] = (('rows', 'cols'), np.zeros(plate_dims))
    Designs['treated_wells'] = (('rows', 'cols'),
                                np.ones(plate_dims, dtype=bool))

    return Designs

