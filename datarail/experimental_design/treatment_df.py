import pandas as pd
import numpy as np


def make_treatment_dataframe(drug_treatment_dict, args,
                             combo_pairs=[], combo_doses=[]):
    """ Function that returns a long table Dataframe for
    drug treatments with n_columns  = len(drugs) and n_rows = n_wells

    Parameters
    ----------
    drug_treatment_dict: dict
             dictionary of drugs as keys and the corresponding doses as values
    args:
       default parameters for plate_dims, stock_concentration
    combo_pairs: list of tuples
              list of drug combinations to be used in the experiment
    combo_doses: dict
             dictionary of drugs as keys and the corresponding doses used in
             combination treatment

    Returns
    -------
    treatment_df: pandas dataframe
             long table dataframe where columns are drugs and rows are wells
    """

    n_wells = np.dot(args.plate_dims[0], args.plate_dims[1])
    total_treatments = len(drug_treatment_dict.keys())
    all_treatments = np.zeros([total_treatments, n_wells])
    count = 0
    drug_num = {drug: i for i, drug in enumerate(drug_treatment_dict.keys())}
    for drug in drug_treatment_dict.keys():
        n_treatments = len(drug_treatment_dict[drug])
        all_treatments[drug_num[drug],
                       count:count+n_treatments] = drug_treatment_dict[drug]
        count += n_treatments
    for pair in combo_pairs:
        n_treatments = len(combo_doses[pair[1]])
        for i in range(len(combo_doses[pair[0]])):
            all_treatments[drug_num[pair[0]],
                           count:count+n_treatments] = combo_doses[pair[0]][i]
            all_treatments[drug_num[pair[1]],
                           count:count+n_treatments] = combo_doses[pair[1]]
        count += n_treatments
    treatment_df = pd.DataFrame(all_treatments.T,
                                columns=drug_treatment_dict.keys())
    return treatment_df
