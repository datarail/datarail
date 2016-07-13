import numpy as np


def randomizer(drugs, xray_struct, all_treatments, cntrl_pos, random_seed):

    plate_dims = cntrl_pos.shape
    nwells = cntrl_pos.size
    cntrl_idx = np.nonzero(cntrl_pos.reshape(1, nwells))[1]
    np.random.seed(random_seed)
    idx = np.random.choice(range(nwells),
                           size=nwells,
                           replace=False)
    idx[cntrl_idx] = nwells
    idx_sort = np.argsort(idx)

    for drug in drugs:
        panel = xray_struct[drug].values
        panel = panel.reshape(1, nwells)
        panel[0, idx_sort[:nwells]] = all_treatments[drug].values
        panel = panel.reshape(plate_dims)
        xray_struct[drug].values = panel
    return xray_struct
