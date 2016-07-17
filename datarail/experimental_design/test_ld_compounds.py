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
    Designs = construct_design(cell_lines, drug_treatment_df,
                               num_doses, args, barcode, random_seed=rep)
    plot_drugs(Designs, drugs, rep)
    plot_control_wells(Designs)


with open('LD_20160713_rep_1.pkl', 'rb') as pickle_file:
    rep1 = pickle.load(pickle_file)
with open('LD_20160713_rep_2.pkl', 'rb') as pickle_file:
    rep2 = pickle.load(pickle_file)
with open('LD_20160713_rep_3.pkl', 'rb') as pickle_file:
    rep3 = pickle.load(pickle_file)


def get_mapping(batch, reps, drug_treatment):
    batch_wells = []
    rep1_wells = []
    rep2_wells = []
    rep3_wells = []
    for i, drug in enumerate(batch):
        row = chr(65+i)
        for j, c in enumerate(drug_treatment[drug]):
            batch_wells.append('%s%d' % (row, j+1))
            rep1_wells += get_well_id(reps[0], drug, c)
            rep2_wells += get_well_id(reps[1], drug, c)
            rep3_wells += get_well_id(reps[2], drug, c)
    z = [batch_wells, rep1_wells, rep2_wells, rep3_wells]
    return z


def get_control_mapping(batch_controls, control_dict, z):
    for control in batch_controls:
        for dose in control_dict[control].keys():
            z[0] += control_dict[control][dose]
            z[1] += get_well_id(reps[0], control, float(dose))
            z[2] += get_well_id(reps[1], control, float(dose))
            z[3] += get_well_id(reps[2], control, float(dose))
    return z        


control_dict = {}
control_dict['DMSO'] = {'300': ['A12', 'B12', 'C12', 'D12', 'E12', 'F12']}
control_dict['Paclitaxel'] = {'300': ['A11', 'B11']}
control_dict['GSK2126458'] = {'300': ['A12', 'B12']}
control_dict['Actinomycin D'] = {'300': ['D11', 'E11', 'F11'], '250': ['G11']}
control_dict['Staurosporine'] = {'300': ['D12', 'E12', 'F12', 'G12']}


batch1_drugs = drugs[:8]
batch2_drugs = drugs[8:16]
batch3_drugs = drugs[16:24]
batch4_drugs = drugs[24:]
batch1_controls = positive_controls
batch2_controls = ['DMSO']

reps = [rep1, rep2, rep3]
columns = ['source_wells', 'destination_wells_rep1',
           'destination_wells_rep2', 'destination_wells_rep3']

z = get_mapping(batch1_drugs, reps, drug_treatment)
z = get_control_mapping(batch1_controls, control_dict, z)
df = pd.DataFrame(zip(z[0], z[1], z[2], z[3]), columns=columns)
df.to_csv('mapping_01.csv', index=False)

z = get_mapping(batch2_drugs, reps, drug_treatment)
z = get_control_mapping(batch2_controls, control_dict, z)
df = pd.DataFrame(zip(z[0], z[1], z[2], z[3]), columns=columns)
df.to_csv('mapping_02.csv', index=False)

z = get_mapping(batch3_drugs, reps, drug_treatment)
df = pd.DataFrame(zip(z[0], z[1], z[2], z[3]), columns=columns)                  
df.to_csv('mapping_03.csv', index=False)

z = get_mapping(batch4_drugs, reps, drug_treatment)
df = pd.DataFrame(zip(z[0], z[1], z[2], z[3]), columns=columns)                  
df.to_csv('mapping_04.csv', index=False)


