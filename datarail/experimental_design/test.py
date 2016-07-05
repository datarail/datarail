from construct_design import construct_design
from treatment_df import get_df
from parse_vararg import parse_args
import numpy as np
from plot_panels import plot_drugs, plot_control_wells

args = parse_args()
drugs = ['MK2206', 'Gefitinib', 'GSK1059615',
         'BEZ235', 'Triciribine', 'PP242']
cells = ['AU565', 'HCC1954', 'T47D', 'MCF 10A']
n_controls = 90

drug_treatment = {}
concentrations = np.hstack((1, 10 ** np.linspace(1, 4, 9)))
for drug in drugs:
    drug_treatment[drug] = concentrations
drug_class1 = ['MK2206', 'Gefitinib', 'GSK1059615']
drug_class2 = ['BEZ235', 'Triciribine', 'PP242']
combo_pairs = [(d1, d2) for d1 in drug_class1 for d2 in drug_class2]
combo_doses = {}
for drug in drugs:
    combo_doses[drug] = concentrations[4:7]
drug_treatment_df = get_df(drug_treatment, args,
                           combo_pairs, combo_doses)

Designs = construct_design(drugs, cells, drug_treatment_df,
                           n_controls, args)
plot_drugs(Designs, drugs)
plot_control_wells(Designs)
