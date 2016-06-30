from construct_design import construct_design
from treatment_df import get_df
from parse_vararg import parse_args
import numpy as np
from plot_drugs import plot

args = parse_args()
drugs = ['MK2206', 'Gefitinib', 'GSK1059615',
         'BEZ235', 'Triciribine', 'PP242']
cells = ['AU565', 'HCC1954', 'T47D', 'MCF 10A']
n_controls = 90

drug_treatment = {}
concentrations = np.hstack((1, 10 ** np.linspace(1, 4, 9)))
for drug in drugs:
    drug_treatment[drug] = concentrations
drug_treatment_df = get_df(drug_treatment, args)

Designs = construct_design(drugs, cells, drug_treatment_df,
                           n_controls, args)
plot(Designs, drugs)
