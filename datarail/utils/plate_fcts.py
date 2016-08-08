import pandas as pd
import xarray as xr

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
    index = pd.MultiIndex.from_arrays([df.barcode, df.row, df.column])
    df.index = index
    df = df.drop(['barcode', 'column', 'row'],axis=1)
    return xr.Dataset.from_dataframe(df)
