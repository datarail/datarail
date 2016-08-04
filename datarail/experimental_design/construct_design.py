from get_hmslid import get_hmslid
from control_position import set_control_positions
from randomizer import randomizer
import xarray as xr
import numpy as np
import cPickle as pickle


def construct_design(treatments_df, args, barcode, n_replicates=1,
                     random_seed=True, edge_bias=True):
    """ construction of plate layout

    Parameters
    ----------
    treatment_df: pandas dataframe
          long table dataframe of drugs and corresponding doses
    args:
        default values of plate dimensions, stock concetntrations
    barcode: str
        barocde for plate
    n_replicates: int
        number of replicates in the experiment
    random_seed: int
        seed for random number generator
    edge_bias: boolean value
        True if edge bias to be taken into account while assigning wells

    Returns
    -------
    Designs: list[xarray replicates]
       list of replicates for Design xarray structures
    """

    treatments = treatments_df.keys()
    # cleaning up the treatments_df to only have rows with a treatment
    drug_treatment_df = treatments_df.iloc[
        (treatments_df.values != 0).any(axis=1), :]
    n_treatments = len(drug_treatment_df)
    plate_dims = args.plate_dims
    plate_rows = list(map(chr, range(65, 65+plate_dims[0])))
    plate_cols = range(1, plate_dims[1] + 1)
    Design1 = xr.Dataset({k: (['rows', 'cols'], np.zeros(args.plate_dims))
                          for k in treatments},
                         coords={'rows': plate_rows, 'cols': plate_cols})

    Design1.attrs['Seed'] = args.Seed
    Design1.attrs['well_volume'] = args.well_volume
    Design1.attrs['plate_dims'] = plate_dims
    # check that the length of barcode match the number of replicates
    # Design1.attrs['barcode'] = barcode

    for treatment in treatments:
        Design1[treatment].attrs['DrugName'] = treatment
        Design1[treatment].attrs['Stock_conc'] = args.stock_conc
        Design1[treatment].attrs['HMSLid'] = get_hmslid([treatment])[treatment]
    Design1['Perturbations'] = (('rows', 'cols'), np.zeros(plate_dims))
    Design1['Vehicle'] = (('rows', 'cols'), np.zeros(plate_dims))

    n_wells = plate_dims[0] * plate_dims[1]
    n_controls = n_wells - n_treatments
    (cntrl_pos, treated_pos) = set_control_positions(
        plate_dims, n_controls)  # forced control positions
    Design1['control_wells'] = (('rows', 'cols'), cntrl_pos)
    Design1['treated_wells'] = (('rows', 'cols'), treated_pos)
    Designs = randomizer(treatments, Design1, drug_treatment_df, plate_dims,
                         cntrl_pos, n_replicates, random_seed, edge_bias)

    for i in range(len(Designs)):
        Designs[i].attrs['barcode'] = barcode[i]
        Designs[i]['control_wells'] = (('rows', 'cols'),
                                       reassign_cntrls(Designs[i], treatments))

    for i in range(len(Designs)):
        filename = '%s.pkl' % barcode[i]
        pickle.dump(Designs, open(filename, 'wb'), protocol=-1)
    return Designs


def reassign_cntrls(Design, treatments):

    all_conc = np.array([Design[d].values for d in treatments])
    return np.all(all_conc == 0, axis=0)


def assign_pc(Design, pc_treatments, well_index):
    """ assign postive control wells to Design plate

    Parameters
    ----------
    Design: xarray sturcture

    pc_treatments: dict
         dict of positive control as keys and doses as values
    well_index: list
         list of wells (ids) that are available for the positive control

    Returns
    -------
    Design: xarray structre
        updates xarray structure with positive controls assigned
    """

    pcs = pc_treatments.keys()
    start = 0
    for pc in pcs:
        Design[pc] = (('rows', 'cols'), np.zeros([16, 24]))
        panel = Design[pc].values
        panel = panel.reshape(1, 384)
        panel[0, well_index[
            start:start+len(pc_treatments[pc])]] = pc_treatments[pc]
        panel = panel.reshape([16, 24])
        Design[pc].values = panel
        start += len(pc_treatments[pc])
    return Design    

