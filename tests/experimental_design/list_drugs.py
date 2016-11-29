import pandas as pd
import numpy as np


df = pd.DataFrame.from_dict({
        'Compound_Name': ['Vehicle', 'D_positive', 'D_fingerprint']+['D_%i'%i for i in range(1,31)],
        'Highest_Dose':['', '1 uM', '5 uM'] + 30*['10 uM'],
        'Role':['negative_control','positive_control','fingerprint']+30*['treatment'],
        'num_wells':[12,6,50]+30*[9]})

df.to_csv('../../examplars/INPUT/compound_list.tsv', sep='\t', index=False)
