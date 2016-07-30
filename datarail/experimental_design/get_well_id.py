import numpy as np


def get_wells(Designs, drug, concentration):
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
