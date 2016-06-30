import pandas as pd
import numpy as np


def get_df(drug_tr_dict, args):
    n_wells = np.dot(args.plate_dims[0], args.plate_dims[1])
    all_treatments = np.zeros([len(drug_tr_dict.keys()), n_wells])
    count = 0
    for i, drug in enumerate(drug_tr_dict.keys()):
        n_treatments = len(drug_tr_dict[drug])
        all_treatments[i, count:count+n_treatments] = drug_tr_dict[drug]
        count += n_treatments

    df = pd.DataFrame(all_treatments.T, columns=drug_tr_dict.keys())
    return df
