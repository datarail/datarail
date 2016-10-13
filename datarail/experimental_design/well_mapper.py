import numpy as np


def get_well_positions(Designs, drug, concentration):
    """ Function that returns the well in a Design xarray
    given the drug name and concentration

    Parameters
    ----------
    Design: xarray structure

    drug: str
      name of drug being queried

    concentration: int
      dose value being queried

    Returns
    -------
    wells: list
       list of wells on the plate
    """
    try:
        arr = Designs[drug].values
    except KeyError:
        print "queried drug not in plate"
    ind = np.where(arr == concentration)
    wells = []
    if len(ind[0]) > 0:
        for i in range(len(ind[0])):
            pos = '%s%s' % (chr(65+ind[0][i]), ind[1][i]+1)
            wells.append(pos)
    return wells


def get_well_index(well, plate_dims):
    """Function that maps well cooridnate to an numerical id
    between 0 and number of wells on the plate

    Parameter
    ---------
    well: str
      well coordinate on plate eg: 'B10'
    plate_dims: np array
       dimensions of plate

    Returns
    -------
    index: int
       numerical id of the well
    """

    row_num = ord(well[0]) - 65
    col_num = int(well[1:]) - 1
    index = row_num*plate_dims[1] + col_num
    return index


def get_inner_untreated_wells(Design, drugs):
    """ function that checks for and returns list of wells
    available in the inner wells

    Parameters
    ----------
    Design: xarray structure
    drugs: list
        list of drugs

    Returns
    -------
    pos_wells: list
         list of well names
    """

    untreated_wells = np.ones([16, 24], dtype=bool)
    treatments = drugs + ['DMSO']
    for i, tr in enumerate(treatments):
        nparray = Design[tr].values
        pos = np.nonzero(nparray)
        for l in range(len(pos[0])):
            untreated_wells[pos[0][l], pos[1][l]] = False
    untreated_wells[0, :] = False
    untreated_wells[-1, :] = False
    untreated_wells[:, 0] = False
    untreated_wells[:, -1] = False
    pos = np.nonzero(untreated_wells)
    rows = [chr(65+i) for i in pos[0]]
    cols = [str(j+1) for j in pos[1]]
    pos_wells = [i+j for i, j in zip(rows, cols)]
    return pos_wells


def get_well_name(well_id, plate_dims):
    """Function that maps well numerical id to name (label)

    Parameter
    ---------
    well_id: int
      well id on plate eg: 234
    plate_dims: np array
       dimensions of plate

    Returns
    -------
    index: str
       coordinate label of the well, eg: 'B10'
    """

    row_num, col_num = divmod(well_id, plate_dims[1])
    row_label = chr(65 + row_num)
    col_label = col_num + 1
    well_name = "%s%02d" % (row_label, col_label)
    return well_name
