from get_hmslid import get_hmslid
from control_position import control_positions
from randomizer import randomizer
import xarray as xr
import numpy as np
import cPickle as pickle


def construct_design(cell_lines, treatments_df,
                     num_doses, args, barcode, random_seed=False):

    treatments = treatments_df.keys()
    plate_dims = args.plate_dims
    plate_rows = list(map(chr, range(65, 65+plate_dims[0])))
    plate_cols = range(1, plate_dims[1] + 1)
    Designs = xr.Dataset({k: (['rows', 'cols'], np.zeros(args.plate_dims))
                          for k in treatments},
                         coords={'rows': plate_rows, 'cols': plate_cols})

    Designs.attrs['Seed'] = args.Seed
    Designs.attrs['well_volume'] = args.well_volume
    Designs.attrs['plate_dims'] = plate_dims
    Designs.attrs['barcode'] = barcode

    for treatment in treatments:
        Designs[treatment].attrs['DrugName'] = treatment
        Designs[treatment].attrs['Stock_conc'] = args.stock_conc
        Designs[treatment].attrs['HMSLid'] = get_hmslid([treatment])[treatment]
    Designs['Perturbations'] = (('rows', 'cols'), np.zeros(plate_dims))
    Designs['Vehicle'] = (('rows', 'cols'), np.zeros(plate_dims))
    Designs['treated_wells'] = (('rows', 'cols'),
                                np.ones(plate_dims, dtype=bool))

    n_wells = plate_dims[0] * plate_dims[1]
    n_controls = n_wells - (len(treatments) * num_doses)
    cntrl_pos = control_positions(plate_dims, n_controls)
    Designs['control_wells'] = (('rows', 'cols'), cntrl_pos)
    Designs = randomizer(treatments, Designs, treatments_df,
                         cntrl_pos, random_seed)
    filename = '%s.pkl' % barcode
    pickle.dump(Designs, open(filename, 'wb'), protocol=-1)
    return Designs
