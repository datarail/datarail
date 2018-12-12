import numpy as np


def get_well_index(well, plate_dims):
    """Function that maps well cooridnate to an numerical id
    between 0 and number of wells on the plate

    Parameters
    ----------
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


def get_well_name(well_id, plate_dims):
    """Function that maps well numerical id to name (label)

    Parameters
    ----------
    well_id : int
      well id on plate eg: 234
    plate_dims : np array
       dimensions of plate

    Returns
    -------
    index : str
       coordinate label of the well, eg: 'B10'
    """

    row_num, col_num = divmod(well_id, plate_dims[1])
    row_label = chr(65 + row_num)
    col_label = col_num + 1
    well_name = "%s%02d" % (row_label, col_label)
    return well_name
