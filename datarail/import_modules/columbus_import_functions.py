import pandas as pd
import numpy as np
import datetime
import re
import os
import datarail.utils.plate_fcts as pltfct



def Columbus_processing(filename, fields,
                        cell_count=None,
                        outputfile=None):
    """ import and clean up the output of Columbus scanner

    Arguments:
    filename -- name of the file containing the data (should be a tab-separated file)
    fields -- list of pairs of names for data to import:
                 column name in the original Columbus file
                 column name in the output file
              e.g.: fields = (('cell_count__total', 'Nuclei-Hoechst Selected - Number of Objects'),
                              ('cell_count__dead', 'Nuclei-LDRpos - Number of Objects'))
    cell_count -- define the arithmetic operation to get the number of viable cell if not defined
                     directly as a field (e.g. 'cell_count__total - cell_count__dead')
    outputfile -- name of the file to save the processed data (default = None, nothing saved)
    """

    df = _Read_Columbus(filename)

    df = _Correct_fields(df)

    df = _Define_count_fields(df, fields, cell_count)

    if outputfile:
        df.to_csv(outputfile, index=False, sep='\t')
        print 'File \n\t%s \nprocessed and saved in \n\t%s' % (filename, outputfile)

    return df


##############################################################

def _Read_Columbus(filename):

    dfin = pd.read_csv(filename, delimiter='\t')

    barcode = [dfin.Result[i].split()[0] for i in range(0,len(dfin))]
    date = [datetime.datetime.strptime(dfin.Result[i].split()[2]+' '+dfin.Result[i].split()[3],
                                       "%Y-%m-%d %H:%M:%S") for i in range(0,len(dfin))]

    dfout = pd.concat([pd.Series(barcode, name='barcode'), pd.Series(date, name='date'), dfin], axis = 1)
    dfout['well'] = [pltfct.well_as_2digit(w) for w in dfout.Well]

    dfout.drop(['Result', 'URL', 'Well'], 1, inplace=True)

    return dfout


def _Correct_fields(dfin):

    dfout = dfin.copy()

    f_cnt = np.bincount(dfin['Number of Analyzed Fields'].astype(int))
    defaultNfield = np.max([i for (i,c) in enumerate(f_cnt) if c==max(f_cnt)])
    print 'Default number of fields: %i ; %i wells with missing field(s)' % \
        (defaultNfield, np.sum(f_cnt[range(0,len(f_cnt)) != defaultNfield]))

    for f in _get_count_fields(dfin):
        dfout[f] = np.ceil((dfin[f]*defaultNfield/dfin['Number of Analyzed Fields']))

    dfout.drop(['Number of Analyzed Fields'], 1, inplace=True)

    return dfout


def _Calculate_time(dfin, time0):
    #### for the timecourse data --- not implemented yet ---
    dfout = dfin.copy()

    return dfout

def _Define_count_fields(dfin,
                         fields,
                         cell_count=None,
                         addfields=[]):

    dfout = dfin.copy()
    fs = _get_count_fields(dfin)

    for af in fields:
        if af[0] not in fs:
            print 'Missing column: ' + af[0]
            raise ValueError('Missing column in dataframe')
        dfout = pd.concat([ dfout, pd.DataFrame(dfin[af[0]].values,
                                                columns=[af[1]], copy=True) ],
                          axis=1)

    if cell_count:
        if 'cell_count' in [fc[1] for fc in fields]:
            print 'cell_count column defined twice'
            raise ValueError('cell_count column defined twice')

        # get the cell count from the fields defined
        cc_f = re.split('\s?[+-/*]\s?', cell_count)
        if np.any([f not in [fc[1] for fc in fields] for f in cc_f]):
            print 'Missing column(s): ' + \
                ' - '.join([f for f in cc_f if f not in [fc[1] for fc in fields]])
            raise ValueError('Missing column in dataframe')

        cc_eval = cell_count
        for f in cc_f:
            cc_eval = re.sub(f, 'dfout.' + f , cc_eval)
        # need to add '.values' in order to cover cases of arithmetic operation AND single value
        dfout = pd.concat([dfout,
                           pd.DataFrame(eval(cc_eval).values, columns=['cell_count'])],
                          axis=1)
        dfout.cell_count[dfout.cell_count<0] = 0


    dfout.drop(fs, 1, inplace=True)


    return dfout



################################## helper functions ###############################

def _get_count_fields(dfin):

    numdf = dfin.copy()
    numdf = numdf.select_dtypes(include=['int16', 'int32', 'int64', 'float16', 'float32', 'float64'])
    for f in ['Row', 'Column', 'Number of Analyzed Fields']:
        if f in numdf.columns:
            numdf.drop(f, 1, inplace=True)

    return [f for f in numdf.columns]
