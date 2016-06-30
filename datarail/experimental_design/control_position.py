import numpy as np


def control_positions(plate_dims, control_count):
    min_inner_positions = 6
    cntrl_pos = np.zeros(plate_dims,
                         dtype=bool)

    assert control_count >= 12,\
        "For placing controls on the edge, atleast 12 controls are required"

    row_pos1 = [round(i-1) for i
                in np.linspace(1, plate_dims[0], 4)]
    row_pos2 = [round(i-1) for i
                in np.linspace(1, plate_dims[0], 3)]
    col_pos = [round(i-1) for i
               in np.linspace(1, plate_dims[1], 4)]

    if control_count < 14:
        print ('Only 6 controls on the edge,'
               'would be better with atleast 14 in controls in total')
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

    if control_count >= (outer_edge + inner_edge + min_inner_positions):
        cntrl_pos[[0, 1, -2, -1], :] = True
        cntrl_pos[:, [0, 1, -2, -1]] = True
    elif control_count >= (outer_edge + min_inner_positions):
        cntrl_pos[[0, -1], :] = True
        cntrl_pos[:, [0, -1]] = True
    elif control_count >= 20:
        cntrl_pos[[0, -1], 0] = True
        cntrl_pos[[0, -1], -1] = True
    return cntrl_pos


def get_boundary_cell_count(plate_dims):
    boundary_cell_count = 2 * (plate_dims[0] + plate_dims[1] - 2)
    return boundary_cell_count
