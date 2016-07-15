

def well_as_1digit(well):
    if well[1]=='0':
        well = well[0]+well[2]
    return well

def well_as_2digit(well):
    if len(well)==2:
        well = well[0]+'0'+well[1]
    return well

def well_as_row_col(well):
    return {'row':ord(well[0])-64,'column':int(well[1:])}
