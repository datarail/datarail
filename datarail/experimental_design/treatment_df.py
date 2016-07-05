import pandas as pd
import numpy as np


def get_df(drug_tr_dict, args, combo_pairs=[], combo_doses=[]):
    n_wells = np.dot(args.plate_dims[0], args.plate_dims[1])
    total_treatments = len(drug_tr_dict.keys())
    all_treatments = np.zeros([total_treatments, n_wells])
    count = 0
    drug_num = {drug: i for i, drug in enumerate(drug_tr_dict.keys())}
    for drug in drug_tr_dict.keys():
        n_treatments = len(drug_tr_dict[drug])
        all_treatments[drug_num[drug],
                       count:count+n_treatments] = drug_tr_dict[drug]
        count += n_treatments
    for pair in combo_pairs:
        n_treatments = len(combo_doses[pair[1]])
        for i in range(len(combo_doses[pair[0]])):
            all_treatments[drug_num[pair[0]],
                           count:count+n_treatments] = combo_doses[pair[0]][i]
            all_treatments[drug_num[pair[1]],
                           count:count+n_treatments] = combo_doses[pair[1]]
        count += n_treatments
    df = pd.DataFrame(all_treatments.T, columns=drug_tr_dict.keys())
    return df

