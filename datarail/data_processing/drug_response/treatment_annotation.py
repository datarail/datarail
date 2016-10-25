import pandas as pd
import numpy as np
import re
from datarail.import_modules.columbus_import_functions import Columbus_processing
import datarail.utils.plate_fcts as pltfct
import warnings


def add_plate_info(dfdata, dfplate):
    ''' add the the metadata from the plate info file'''
    # need to implement tests on the imputs

    dfout = pd.merge(dfdata, dfplate, on='barcode')

    return dfout


def add_treatments(dfdata, trtfolder):
    ''' add the the data about the treatment'''

    trt_files = list(set(dfdata.treatment_file) - set(['-', '']))
    # add the extension if not present (.tsv by default)
    trt_filenames = [f+'.tsv' if (re.search('\.\w\w\w$', f) is None) else f for f in trt_files]

    trt = [pd.read_csv(trtfolder + '/' + f, sep='\t') for f in trt_filenames]
    for i in range(len(trt)):
        trt[i]['well'] = [pltfct.well_as_2digit(w) for w in trt[i]['well']]

    dfout = pd.merge(dfdata.loc[np.any([dfdata.treatment_file == '',dfdata.treatment_file =='-'], axis=0),:],
                     _untreated_plate(), on='well')
    for (i,tf) in enumerate(trt_files):
        dfout = dfout.append(pd.merge(dfdata.loc[dfdata.treatment_file == tf, :],
                                      trt[i], on='well'))

    ### error message if some wells have been dropped
    if len(dfout) != len(dfdata):
        warnings.warn('wells have been dropped, output length is different than input length', stacklevel=2)

    return dfout


def average_replicates(dfdata, keys=None):
    if keys is None:
        keys = set(['agent', 'concentration', 'cell_line', 'treatment_duration', 'role']).intersection(dfdata.columns)
    elif not all([keys[i] in dfdata.columns for i in range(len(keys))]):
        raise Exception('Not all keys are a column of data table')
    else:
        keys = set(keys)

    df_ = dfdata.copy()
    if 'concentration' in keys:
        # rounding of concentration to avoid issues due to precision
        df_.concentration = (10**np.round(np.log10(df_.concentration),3))

    # find the columns that will be averaged
    numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
    numcol = set(df_.select_dtypes(include = numerics).columns) -keys -set(['Column', 'Row'])

    print 'Columns to average: "%s"' % '" "'.join(numcol)

    # find the annotation columns
    c_annot = df_.loc[:,list(set(df_.columns) -numcol -keys -set(['Column', 'Row']))].columns.values

    c_l = [c for c in c_annot if len(df_.loc[:,list(keys) + [c]].drop_duplicates())==\
           len(df_.loc[:,list(keys)].drop_duplicates())]

    if len(c_l)>0:
        print 'Columns added as annotations: "%s"' % '" "'.join(c_l)
    if len(c_annot)>len(c_l):
        print('\n-->Following columns are discarded:\n "%s" \n\t(set as key if necessary)' %
                      '" "'.join(list(set(c_annot) - set(c_l))))


    dfout = df_.loc[:,list(keys) + list(numcol)].groupby(list(keys), as_index=False).mean()
    dfout = pd.merge(dfout, df_.loc[:,list(keys) + c_l].drop_duplicates(), on=list(keys))

    return dfout.loc[:,list(keys) + c_l + list(numcol)]


def load_example():
    dfdata = Columbus_processing('../../drug_response_data/OUTPUT/Example1_Columbus_output.tsv')

    dfplate = pd.read_csv('../../drug_response_data/OUTPUT/Example1_plate_info.tsv', sep='\t')
    trtfolder = '../../drug_response_data/OUTPUT/'

    return {'data':dfdata, 'plates':dfplate, 'folder':trtfolder}

def _untreated_plate():
    well, agent, concentration, role = [[] for _ in range(4)]
    for ir in range(1,17):
        for ic in range(1,25):
            well += [chr(64+ir)+str(ic).zfill(2)]
            agent += ['-']
            concentration += [0]
            role += ['untreated']
    return pd.DataFrame.from_dict(dict({
        'well':well,
        'agent':agent,
        'concentration':concentration,
        'role':role}))
