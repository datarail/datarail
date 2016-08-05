import pandas as pd
import numpy as np
import datetime
import re
import os
from datarail.utils.plate_fcts import *


def Read_Columbus(filename):

    dfin = pd.read_csv(filename, delimiter='\t')

    barcode = [dfin.Result[i].split()[0] for i in range(0,len(dfin))]
    date = [datetime.datetime.strptime(dfin.Result[i].split()[2]+' '+dfin.Result[i].split()[3],
                                       "%Y-%m-%d %H:%M:%S") for i in range(0,len(dfin))]

    dfout = pd.concat([pd.Series(barcode, name='barcode'), pd.Series(date, name='date'), dfin], axis = 1)
    dfout['well'] = [well_as_2digit(w) for w in dfout.Well]

    dfout.drop(['Result', 'URL', 'Well'], 1, inplace=True)


    return dfout


def Correct_fields(dfin):

    dfout = dfin.copy()

    f_cnt = np.bincount(dfin['Number of Analyzed Fields'].astype(int))
    defaultNfield = np.max([i for (i,c) in enumerate(f_cnt) if c==max(f_cnt)])
    print 'Default number of fields: %i ; %i wells with missing field(s)' % \
        (defaultNfield, np.sum(f_cnt[range(0,len(f_cnt)) != defaultNfield]))

    for f in get_count_fields(dfin):
        dfout[f] = np.ceil((dfin[f]*defaultNfield/dfin['Number of Analyzed Fields']))

    dfout.drop(['Number of Analyzed Fields'], 1, inplace=True)

    return dfout


def Calculate_time(dfin, time0):
    #### for the timecourse data --- not implemented yet ---
    dfout = dfin.copy()

    return dfout

def Define_count_fields(dfin,
                        cell_count='Hoechst_pos-Hoechst_LDR_pos',
                        addfields=[]):

    fs = get_count_fields(dfin)

    # get the cell count from the fields from Columbus
    cc_f = re.split('[+-/*]', cell_count)
    if np.any([f not in fs for f in cc_f]):
        print 'Missing column(s): ' + \
            ' - '.join([f for f in cc_f if f not in fs])
        raise ValueError('Missing column in dataframe')

    cc_eval = cell_count
    for f in cc_f:
        cc_eval = re.sub(f, 'dfin.' + f, cc_eval)

    # need to add '.values' in order to cover cases of arithmetic operation AND single value
    cc = pd.DataFrame(eval(cc_eval).values, columns=['cell_count'])

    for af in addfields:
        a_f = re.split('[+-/*]', af[1])

        if np.any([f not in fs for f in a_f]):
            print 'Missing column(s): ' + \
                ' - '.join([f for f in cc_f if f not in fs])
            raise ValueError('Missing column in dataframe')

        a_eval = af[1]
        for f in a_f:
            a_eval = re.sub(f, 'dfin.' + f, a_eval)
        # need to add '.values' in order to cover cases of arithmetic operation AND single value
        cc = pd.concat([cc, pd.DataFrame(eval(a_eval).values, columns=[af[0]])],
                       axis=1)


    dfout = pd.concat([dfin, cc], axis=1)
    dfout.drop(fs, 1, inplace=True)

    return dfout


def Columbus_processing(filename, outputfile=''):

    df = Read_Columbus(filename)

    df = Correct_fields(df)

    df = Define_count_fields(df,
                        cell_count='Hoechst_pos-Hoechst_LDR_pos',
                        addfields=(('cell_count__total', 'Hoechst_pos'),
                                   ('corpse_count','LDR_pos_Hoechst_neg'),
                                   ('cell_count__dead', 'Hoechst_LDR_pos')))

    if len(outputfile)!=0:
        df.to_csv(outputfile, index=False, sep='\t')
        print 'File \n\t%s \nprocessed and saved in \n\t%s' % (filename, outputfile)

    return df


################################## helper functions ###############################

def get_count_fields(dfin):

    numdf = dfin.copy()
    numdf = numdf.select_dtypes(include=['int16', 'int32', 'int64', 'float16', 'float32', 'float64'])
    for f in ['Row', 'Column', 'Number of Analyzed Fields']:
        if f in numdf.columns:
            numdf.drop(f, 1, inplace=True)

    return [f for f in numdf.columns if np.all(numdf[f]==np.ceil(numdf[f]))]
