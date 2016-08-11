# from datarail.databases import lincs_client
import xarray as xr
import numpy as np
# import cPickle as pickle
from process_assay import make_treatment_dataframe
import edge_barcode
import well_mapper


def make_layout(treatments_dict, barcode_prefix, encode_barcode=True,
                plate_dims=[16, 24], nreps=1, randomize=True,
                biased_randomization=True,
                combo_pairs=[], combo_doses=[], combo_k=1,
                Seed=50, well_volume=10, stock_conc='50uM'):
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
    treatments_df = make_treatment_dataframe(treatments_dict, plate_dims,
                                             combo_pairs, combo_doses)
    treatment_df_wells = assign_well_index(treatments_df, plate_dims, nreps,
                                           randomize, biased_randomization)
    
    barcodes = [barcode_prefix + chr(65+i) for i in range(nreps)]
    Designs = []
    for rep in range(1, nreps+1):
        tr_panel, conc_panel = make_arrays(treatment_df_wells, plate_dims,
                                           combo_k, rep)
        if encode_barcode:
            bc_treatments = treatments_dict[3]
            tr_panel, conc_panel = assign_bc(bc_treatments, barcodes, rep,
                                             tr_panel, conc_panel,
                                             combo_k, plate_dims)
        design = make_xr(conc_panel, tr_panel, plate_dims, combo_k)
        Designs.append(design)
    return Designs


def reassign_cntrls(Design, treatments):

    all_conc = np.array([Design[d].values for d in treatments])
    return np.all(all_conc == 0, axis=0)


def assign_bc(bc_treatments, barcodes, rep, tr_panel, conc_panel,
              combo_k, plate_dims):
    rep_ind = rep - 1
    bc_wells = edge_barcode.encode_barcode(barcodes[rep_ind])
    well_index = [well_mapper.get_well_index(well, plate_dims=[16, 24])
                  for well in bc_wells]
    print well_index
    conc_list, bc_list = [], []
    for k, v in bc_treatments.iteritems():
        bc_list += [k]*len(v)
        conc_list += v
    n_wells = plate_dims[0] * plate_dims[1]
    tr_panel = tr_panel.reshape([1, n_wells, combo_k])
    conc_panel = conc_panel.reshape([1, n_wells, combo_k])

    for id, well_id in enumerate(well_index):
        tr_panel[0, well_id, 0] = bc_list[id]
        conc_panel[0, well_id, 0] = conc_list[id]
    conc_panel = conc_panel.reshape([plate_dims[0], plate_dims[1], combo_k])
    tr_panel = tr_panel.reshape([plate_dims[0], plate_dims[1], combo_k])
    return tr_panel, conc_panel


def edge_bias_randomizer(n_treatments, cntrl_idx,
                         plate_dims, n_replicate):

    np.random.seed(1)
    n_edge = 5
    # split the wells based on distance from edge
    well_groups = edge_wells(plate_dims, n_edge)
    # remove fixed controls
    well_groups = [list(set(w)-set(cntrl_idx)) for w in well_groups]
    # number of wells available
    n_wells = sum([len(w) for w in well_groups])
    n_cntrls = [int(np.floor((n_wells-n_treatments)*(len(w)*1./n_wells)))
                for w in well_groups]
    # assign the number of controls foe each layer
    n_cntrls[-1] = n_wells - n_treatments - sum(n_cntrls[:-1])
    n_trtwells = [len(w) - n_cntrls[i] for i, w in enumerate(well_groups)]

    # count how often a condition is at the distance i from the edge
    edge_cnt = np.zeros([n_treatments, n_edge])
    trt_positions = -np.ones([n_replicate, n_treatments])

    for i_rep in range(n_replicate):
        min_cnt = edge_cnt.min(0)
        # bypass the condition for the inner most set of wells
        min_cnt[-1] = n_replicate

        for j in (i for i in range(n_edge) if n_trtwells[i] > 0):
            # going through each layer and finding treatment to assign

            # picking the positions (left ones will be the controls)
            pos = np.random.choice(range(len(well_groups[j])),
                                   n_trtwells[j], replace=False)
            pos = [w for i, w in enumerate(well_groups[j]) if i in pos]

            # picking the treatments
            candidate_trts = np.nonzero(np.logical_and(
                edge_cnt[:, j] <= min_cnt[j],
                trt_positions[i_rep, :] == -1))[0]

            if len(candidate_trts) >= n_trtwells[j]:
                idx = np.random.choice(range(len(candidate_trts)),
                                       n_trtwells[j], replace=False)
            else:
                idx = candidate_trts
                candidate_trts = np.nonzero(np.logical_and(
                    edge_cnt[:, j] <= min_cnt[j]+1,
                    trt_positions[i_rep, :] == -1)[0])
                # assuming that the number of treated wells will
                # be filled up on the second round
                # -- may need an iterative approach
                idx += np.random.choice(range(len(candidate_trts)),
                                        n_trtwells[j]-len(idx), replace=False)

            # print candidate_trts
            # print (i_rep, j, len(candidate_trts), len(pos), len(idx))
            # print np.vstack((np.array(pos),idx, candidate_trts[idx])).T

            edge_cnt[candidate_trts[idx], j] += 1
            trt_positions[i_rep, candidate_trts[idx]] = pos
            # print trt_positions[i_rep,:]
            # print (sum(edge_cnt[:, j]),
            #       sum(trt_positions[i_rep, :] == -1),
            #       sum(trt_positions[i_rep, :] >= 0))

        assert(np.all(trt_positions[i_rep, :] >= 0))

    return trt_positions


def unbiased_randomizer(n_treatments, cntrl_idx,
                        plate_dims, n_replicate=1):
    np.random.seed(1)
    trt_positions = -np.ones([n_replicate, n_treatments])
    wells = list(set(range(plate_dims[0]*plate_dims[1])) - set(cntrl_idx))
    for i_rep in range(n_replicate):
        pos = np.random.choice(wells, size=n_treatments, replace=False)
        trt_positions[i_rep, :] = pos
    return trt_positions


def edge_wells(plate_dims, n_edge):
    n_wells = plate_dims[0]*plate_dims[1]
    well_dist = np.array([[np.min([i, plate_dims[1]-i-1, j,
                                   plate_dims[0]-j-1])
                           for i in range(plate_dims[1])]
                          for j in range(plate_dims[0])])
    well_groups = [[j for j, w in enumerate(
        well_dist.reshape(n_wells, 1)) if w[0] == i]
                   for i in range(n_edge-1)]
    well_groups.append([i for i in range(n_wells)
                        if i not in sum(well_groups, [])])
    return well_groups


def set_control_positions(plate_dims, control_count):
    """ Number and position of control wells assigned based on
    plate dimension and number of available wells for control

    Parameters
    ---------
    plate_dims: array
         dimension of plate
    control_count: int
         number of wells available for controls

    Returns
    -------
    control_pos: boolean array
          array of dimensions mathing plate dims which
         are True for wells assigned as control wells
    """

    min_inner_positions = 6
    cntrl_pos = np.zeros(plate_dims,
                         dtype=bool)
    treated_pos = np.ones(plate_dims,
                          dtype=bool)

    assert control_count >= 12,\
        "For placing controls on the edge, at least 12 controls are required"

    row_pos1 = [int(i-1) for i
                in np.linspace(1, plate_dims[0], 4)]
    row_pos2 = [int(i-1) for i
                in np.linspace(1, plate_dims[0], 3)]
    col_pos = [int(i-1) for i
               in np.linspace(1, plate_dims[1], 4)]

    # fixed positions for the controls
    if control_count < 14:
        print ('Only 6 controls on the edge,'
               'would be better with at least 14 in controls in total')
        cntrl_pos[row_pos1, col_pos[1]] = True
        cntrl_pos[row_pos1, col_pos[2]] = True
        cntrl_pos[row_pos2[1], col_pos] = True
    else:
        cntrl_pos[row_pos1 + [row_pos2[1]], col_pos[1]] = True
        cntrl_pos[row_pos1 + [row_pos2[1]], col_pos[2]] = True
        cntrl_pos[row_pos1[1], col_pos] = True
        cntrl_pos[row_pos1[2], col_pos] = True

    outer_edge = get_boundary_cell_count(plate_dims)
    inner_dims = [i-2 for i in plate_dims]
    inner_edge = get_boundary_cell_count(inner_dims)

    # remove the outer wells if enough controls
    if control_count >= (outer_edge + inner_edge + min_inner_positions):
        # first and second outer layers
        cntrl_pos[[0, 1, -2, -1], :] = True
        cntrl_pos[:, [0, 1, -2, -1]] = True
        treated_pos[[0, 1, -2, -1], :] = False
        treated_pos[:, [0, 1, -2, -1]] = False

    elif control_count >= (outer_edge + min_inner_positions):
        # outer most layer
        cntrl_pos[[0, -1], :] = True
        cntrl_pos[:, [0, -1]] = True
        treated_pos[[0, -1], :] = False
        treated_pos[:, [0, -1]] = False
    elif control_count >= 20:
        # only the corners
        cntrl_pos[[0, -1], 0] = True
        cntrl_pos[[0, -1], -1] = True
        treated_pos[[0, -1], 0] = False
        treated_pos[[0, -1], -1] = False

    return (cntrl_pos, treated_pos)


def get_boundary_cell_count(plate_dims):
    """ get number of wells in outer or inner edges

    Parameter
    ---------
    plate_dims: array
         dimensions of plate

    Returns
    -------
    boundary_cell_count: int
           number of wells in the edges
    """
    boundary_cell_count = 2 * (plate_dims[0] + plate_dims[1] - 2)
    return boundary_cell_count


def assign_well_index(treatment_df, plate_dims, n_replicates,
                      randomize=False, biased_randomization=False):
    """ Given treatment dataframe, return additional column
        with randmized well id

    Parameter
    ---------
    treatment_df: pandas dataframe
         long table of treatments as columns and doses as rows

    plate_dims: list
         physical dimensions of the plate

    n_replicates: int
         number of replicates

    randomize: boolean
         True if allotment of drug-dose to wells is randomized

    biased_randomization: boolean
         True if randomization is biased in order to avoid assigning
         a given drug-dose combination to the same well across replicates

    Returns
    -------
    dfw_wells: pandas dataframe
          treatment_df appended with the addtional column of well index

    """

    dfw_wells = treatment_df.copy()
    n_wells = plate_dims[0] * plate_dims[1]
    n_treatments = len(dfw_wells)
    n_controls = n_wells - n_treatments
    cntrl_pos, _ = set_control_positions(plate_dims, n_controls)
    cntrl_idx = np.nonzero(cntrl_pos.reshape(1, n_wells))[1]
    if randomize:
        if biased_randomization:
            # randomization with contol of edge locations
            # (keeping fixed controls at their place)
            trt_positions = edge_bias_randomizer(n_treatments, cntrl_idx,
                                                 plate_dims, n_replicates)
        else:
            # unbiased randomization (keeping fixed controls at their place)
            trt_positions = unbiased_randomizer(n_treatments, cntrl_idx,
                                                plate_dims, n_replicates)
    else:
        # no randomization: listing conditions skipping fixed controls
        wells = list(set(range(plate_dims[0]*plate_dims[1])) - set(cntrl_idx))
        trt_positions = [np.array(wells[:n_treatments])
                         for i in range(n_replicates)]

    for i in range(n_replicates):
        col_name = 'rep%d_well_index' % (i+1)
        dfw_wells[col_name] = trt_positions[i, :]
    return dfw_wells


def make_arrays(treatment_df_wells, plate_dims, combo_k, rep):
    """ makes 3d array of concentrations, treatments that will be used
    as input to construct the xarray data structure

    Parameters
    ----------
    treatment_df_wells: pandas dataframe
        long table of treatments

    plate_dims: list
       physical plate dimensions

    combo_k: int
       the max number of treatments in a single well

    rep: int
      the replicate number

    Returns
    -------
    tr_panel: 3-dimensional array
       array of treaments
    conc_panel, 3-dimension array
       array of treatment concentrations
    """

    df = treatment_df_wells.copy()
    cols = df.columns.tolist()
    treatments = [t for t in cols if 'well_index' not in t]
    query_column = 'rep%d_well_index' % rep
    indeces = df[query_column].tolist()
    n_wells = plate_dims[0] * plate_dims[1]
    conc_panel = np.zeros([1, n_wells, combo_k])
    tr_panel = np.zeros([1, n_wells, combo_k], dtype='|S20')

    for ind in indeces:
        df_row = df[df[query_column] == ind]
        tr_list = df_row[df_row[treatments] != 0].dropna(
            axis=1).columns.tolist()
        conc_list = df_row[df_row[treatments] != 0].dropna(
            axis=1).values[0].tolist()
        if len(tr_list) == 1:
            tr_panel[0, ind, 0] = tr_list[0]
            conc_panel[0, ind, 0] = conc_list[0]
        elif len(tr_list) > 1:
            for tr in range(1, len(tr_list)+1):
                print tr
                tr_panel[0, ind, tr] = tr_list[tr-1]
                conc_panel[0, ind, tr] = conc_list[tr-1]
    tr_panel = tr_panel.reshape([16, 24, combo_k])
    conc_panel = conc_panel.reshape([16, 24, combo_k])
    return tr_panel, conc_panel


def make_xr(conc_panel, tr_panel, plate_dims, combo_k):
    plate_rows = list(map(chr, range(65, 65+plate_dims[0])))
    plate_cols = range(1, plate_dims[1] + 1)
    combo_dims = ['single']
    if combo_k > 1:
        combo_dims += ["combo_drug_%d" % d for d in range(1, combo_k)]

    design = xr.Dataset({'concentrations': (['rows', 'cols', 'combos'],
                                            conc_panel)},
                        coords={'treatments': (['rows', 'cols', 'combos'],
                                               tr_panel),
                                'rows': plate_rows, 'cols': plate_cols,
                                'combos': combo_dims})
    return design
