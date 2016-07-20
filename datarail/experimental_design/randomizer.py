import numpy as np
import xarray as xr


def randomizer(drugs, xray_struct, all_treatments,
               plate_dims, cntrl_pos, n_replicates,
               random_seed=1, edge_bias=True):

    n_wells = plate_dims[0]*plate_dims[1]
    n_treatments = len(all_treatments)
    cntrl_idx = np.nonzero(cntrl_pos.reshape(1, n_wells))[1]
    if random_seed:
        np.random.seed(random_seed)
        if edge_bias:
            # randomization with contol of edge locations
            # (keeping fixed controls at their place)
            trt_positions = edge_bias_randomizer(n_treatments, cntrl_idx,
                                                 random_seed, plate_dims,
                                                 n_replicates)
        else:
            # unbiased randomization (keeping fixed controls at their place)
            trt_positions = unbiased_randomizer(n_treatments, cntrl_idx,
                                                random_seed, plate_dims,
                                                n_replicates)
    else:
        # no randomization: listing conditions skipping fixed controls
        wells = list(set(range(plate_dims[0]*plate_dims[1])) - set(cntrl_idx))
        trt_positions = [np.array(wells[:n_treatments])
                         for i in range(n_replicates)]

    xray_structs = [xr.Dataset() for i in range(n_replicates)]
    for i in range(n_replicates):
        xray_structs[i] = xray_struct.copy(deep=True)

        for drug in drugs:
            panel = xray_structs[i][drug].values
            panel = panel.reshape(1, n_wells)
            panel[0, trt_positions[i].astype('int')] = all_treatments[drug].values
            panel = panel.reshape(plate_dims)
            xray_structs[i][drug].values = panel

    return xray_structs


def edge_bias_randomizer(n_treatments, cntrl_idx,
                         random_seed, plate_dims, n_replicate):

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
                        random_seed, plate_dims, n_replicate=1):

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
                           for i in range(plate_dims[1])] for j in range(plate_dims[0])])
    well_groups = [[j for j,w in enumerate(well_dist.reshape(n_wells,1)) if w[0] == i]
                   for i in range(n_edge-1)]
    well_groups.append([i for i in range(n_wells)
                        if i not in sum(well_groups, [])])
    return well_groups
