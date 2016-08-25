import pandas as pd
import xarray as xr
import numpy as np

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

def add_row_col(df):
    return pd.concat([df,
                    pd.DataFrame( [well_as_row_col( df.well[i] )['row']
                                   for i in range(len(df))], columns=['row']),
                    pd.DataFrame( [well_as_row_col( df.well[i] )['column']
                                   for i in range(len(df))], columns=['column'])],
                   axis=1)

def dfplate2xr(df):
    df = add_row_col(df)

    # deduce the size of the plate based on standard format (2^n x 1.5*2^n)
    n = max(np.ceil(np.log2(max(df.column)/1.5)), np.ceil(np.log2(max(df.row))))

    # clean the df to convert to xr
    index = pd.MultiIndex.from_arrays([df.barcode, df.row, df.column])
    df.index = index
    df = df.drop(['barcode', 'column', 'row'],axis=1)

    xray =  xr.Dataset.from_dataframe(df)

    xray.attrs['plate_dims'] = [int(2**n), int(1.5*2**n)]

    return xray

def axis_plate(h, plate_dims):

    h.axes.axes.set_ylim([.5, plate_dims[0]+.5])
    h.axes.axes.set_xlim([.5, plate_dims[1]+.5])

    if plate_dims[0]>10:
        h.axes.set_yticks(range(1, plate_dims[0]+1,2))
        h.axes.set_yticklabels([chr(i) for i in range(ord('A'), ord('A')+plate_dims[0], 2)])
    else:
        h.axes.set_yticks(range(1,plate_dims[0]+1))
        h.axes.set_yticklabels([chr(i) for i in range(ord('A'), ord('A')+plate_dims[0])])

    if plate_dims[1]>15:
        h.axes.set_xticks(range(1,plate_dims[1]+1,2))
    else:
        h.axes.set_xticks(range(1,plate_dims[1]+1))

    h.axes.invert_yaxis()
