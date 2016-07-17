from construct_design import construct_design
from treatment_df import get_df
from parse_vararg import parse_args
import numpy as np
from plot_panels import plot_drugs, plot_control_wells
import pandas as pd
import re
from collections import OrderedDict
import cPickle as pickle
from get_well_id import get_well_id

args = parse_args()

df = pd.read_csv('laura_compound_list.tsv', sep='\t')
compounds = df.Compound_Name.tolist()
drugs = compounds[:-4]
positive_controls = compounds[-4:]

cell_lines = ['SKBR3']
num_doses = 9
num_replicates = 3
plate_name = 'LD_20160713'

drug_treatment = OrderedDict()
for drug in drugs:
    max_dose = df['Highest_Dose'].ix[
        df['Compound_Name'] == drug].values[0]
    max_dose_value = float(re.split('u', max_dose)[0])
    drug_treatment[drug] = max_dose_value * 1e-4 * np.logspace(0, 4, num_doses)
drug_treatment['Paclitaxel'] = np.array([300]*2)
drug_treatment['GSK2126458'] = np.array([300]*2)
drug_treatment['Staurosporine'] = np.array([300]*4)
drug_treatment['Actinomycin D'] = np.array([300, 300, 300, 250])
drug_treatment['DMSO'] = np.array([300]*6)
drug_treatment_df = get_df(drug_treatment, args)

for rep in range(num_replicates):
    rep += 1
    barcode = plate_name + '_rep_%d' % rep
    Designs = construct_design(drugs, cell_lines, drug_treatment_df,
                               num_doses, args, barcode, random_seed=rep)
    plot_drugs(Designs, drugs, rep)
    plot_control_wells(Designs)


with open('LD_20160713_rep_1.pkl', 'rb') as pickle_file:
    rep1 = pickle.load(pickle_file)
with open('LD_20160713_rep_2.pkl', 'rb') as pickle_file:
    rep2 = pickle.load(pickle_file)
with open('LD_20160713_rep_3.pkl', 'rb') as pickle_file:
    rep3 = pickle.load(pickle_file)


batch1_drugs = drugs[:8]
batch2_drugs = drugs[8:16]
batch3_drugs = drugs[16:24]
batch4_drugs = drugs[24:]


batches = [batch1_drugs, batch2_drugs,
           batch3_drugs, batch4_drugs]
reps = [rep1, rep2, rep3]


def get_mapping(batch, reps, drug_treatment):
    batch_wells = []
    rep1_wells = []
    rep2_wells = []
    rep3_wells = []
    for i, drug in enumerate(batch):
        row = chr(65+i)
        for j, c in enumerate(drug_treatment[drug]):
            batch_wells.append('%s%d' % (row, j))
            rep1_wells.append(get_well_id(reps[0], drug, c))
            rep2_wells.append(get_well_id(reps[1], drug, c))
            rep3_wells.append(get_well_id(reps[2], drug, c))
    z = zip(batch_wells, rep1_wells, rep2_wells, rep3_wells)
    df = pd.DataFrame(z, columns=['batch1', 'rep1', 'rep2', 'rep3'])
    return df


for i, batch in enumerate(batches):
    df = get_mapping(batch, reps, drug_treatment)
    filename = 'mapping_0%d.csv' % (i+1)
    df.to_csv(filename, index=False)
