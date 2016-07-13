from construct_design import construct_design
from treatment_df import get_df
from parse_vararg import parse_args
import numpy as np
from plot_panels import plot_drugs, plot_control_wells
import pandas as pd
import re

args = parse_args()

df = pd.read_csv('laura_compound_list.tsv', sep='\t')
compounds = df.Compound_Name.tolist()
drugs = compounds[:-4]
positive_controls = compounds[-4:]

cell_lines = ['SKBR3']
num_doses = 9
num_replicates = 3
plate_name = 'LD_20160713'

drug_treatment = {}
for drug in drugs:
    max_dose = df['Highest_Dose'].ix[
        df['Compound_Name'] == drug].values[0]
    max_dose_value = float(re.split('u', max_dose)[0])
    drug_treatment[drug] = max_dose_value * 1e-4 * np.hstack(
        (1, 10 ** np.linspace(1, 4, num_doses)))
drug_treatment_df = get_df(drug_treatment, args)

for rep in range(num_replicates):
    rep += 1
    barcode = plate_name + '_rep_%d' % rep
    Designs = construct_design(drugs, cell_lines, drug_treatment_df,
                               num_doses, args, barcode, random_seed=rep)
    plot_drugs(Designs, drugs, rep)
    plot_control_wells(Designs)

    
