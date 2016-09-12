import pandas as pd
import numpy as np
import re
from datarail.import_modules.columbus_import_functions import Columbus_processing




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

    dfout = pd.merge(dfdata.loc[np.any([dfdata.treatment_file == '',dfdata.treatment_file =='-'], axis=0),:],_untreated_plate(), on='well')
    for (i,tf) in enumerate(trt_files):
        dfout = dfout.append(pd.merge(dfdata.loc[dfdata.treatment_file == tf, :],
                                      trt[i], on='well'))

    return dfout


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
