import pandas as pd
import numpy as np
from datarail.import_modules.columbus_import.columbus_import_functions import Columbus_processing




def add_plate_info(dfdata, dfplate):
    ''' add the the metadata from the plate info file'''
    # need to implement tests on the imputs

    dfout = pd.merge(dfdata, dfplate, on='barcode')

    return dfout


def add_treatments(dfdata, trtfolder):
    ''' add the the data about the treatment'''

    trt_files = dfdata.treatment_file.unique()
    trt = [pd.read_csv(trtfolder + f, sep='\t') for f in trt_files]

    dfout = pd.DataFrame([])
    for (i,tf) in enumerate(trt_files):
        dfout = dfout.append(pd.merge(dfdata.loc[dfdata.treatment_file == tf, :],
                                      trt[i], on='well'))

    return dfout


def load_example():
    dfdata = Columbus_processing('../../drug_response_data/OUTPUT/Example1_Columbus_output.tsv')

    dfplate = pd.read_csv('../../drug_response_data/OUTPUT/Example1_plate_info.tsv', sep='\t')
    trtfolder = '../../drug_response_data/OUTPUT/'

    return {'data':dfdata, 'plates':dfplate, 'folder':trtfolder}
