def trimmean(vals):

    return np.mean(sorted(vals)[int(len(vals)*.25):int(len(vals)*.75)])


def well_as_1digit(well):
    if well[1]=='0':
        well = well[0]+well[2]
    return well


def well_as_row_col(well):
    return {'row':ord(well[0])-64,'column':int(well[1:])}
