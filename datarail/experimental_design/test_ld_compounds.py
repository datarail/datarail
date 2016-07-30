from construct_design import construct_design
from treatment_df import get_df
from parse_vararg import parse_args
import numpy as np
from plot_panels import plot_drugs, plot_control_wells
import pandas as pd
import re
from collections import OrderedDict
import cPickle as pickle
import get_well_id
import edge_barcode

args = parse_args()

df = pd.read_csv('laura_compound_list.tsv', sep='\t')
compounds = df.Compound_Name.tolist()
drugs = compounds[:-4]
positive_controls = compounds[-4:]

cell_lines = ['SKBR3']
num_doses = 9
num_replicates = 3
barcodes = ['DUBI1_' + chr(65+i) for i in range(num_replicates)]

drug_treatment = OrderedDict()
for drug in drugs:
    max_dose = df['Highest_Dose'].ix[
        df['Compound_Name'] == drug].values[0]
    max_dose_value = float(re.split('u', max_dose)[0])
    drug_treatment[drug] = max_dose_value * 1e-4 * np.logspace(0, 4, num_doses)
drug_treatment['DMSO'] = np.array([300]*12)
drug_treatment_df = get_df(drug_treatment, args)

pc_treatments = OrderedDict()
pc_treatments['pc_Staurosporine'] = np.array([1]*8)
pc_treatments['pc_Actinomycin_D'] = np.array([1]*8)
pc_treatments['pc_Paclitaxel'] = np.array([1]*4)
pc_treatments['pc_GSK2126458'] = np.array([0.1]*4)

Designs = construct_design(drug_treatment_df, args, barcodes, 
                           n_replicates=len(barcodes),
                           random_seed=1, edge_bias=True)

for design in Designs:
    bc_wells = edge_barcode.encode_barcode(design.barcode)
    inner_empty_wells = edge_barcode.get_inner_untreated_wells(design, drugs)
    all_wells = bc_wells + inner_empty_wells
    well_index = [get_well_id.wellid2index(well, plate_dims=[16, 24])
                  for well in all_wells]
    design = edge_barcode.assign_pc(design, pc_treatments, well_index)
    # Assign 'DMSO' to remainder wells
    inner_empty_wells = edge_barcode.get_inner_untreated_wells(design, drugs)
    well_index = [get_well_id.wellid2index(well, plate_dims=[16, 24])
                  for well in inner_empty_wells]
    panel = design['DMSO'].values
    panel = panel.reshape(1, 384)
    panel[0, well_index] = 300
    design['DMSO'].values = panel.reshape([16, 24])

reps = Designs


def get_mapping(batch, reps, drug_treatment):
    batch_wells = []
    rep1_wells = []
    rep2_wells = []
    rep3_wells = []
    for i, drug in enumerate(batch):
        row = chr(65+i)
        for j, c in enumerate(drug_treatment[drug]):
            batch_wells.append('%s%d' % (row, j+1))
            rep1_wells += get_well_id.get_well_id(reps[0], drug, c)
            rep2_wells += get_well_id.get_well_id(reps[1], drug, c)
            rep3_wells += get_well_id.get_well_id(reps[2], drug, c)
    z = [batch_wells, rep1_wells, rep2_wells, rep3_wells]
    return z


def get_control_mapping(batch_controls, control_dict, z):
    for control in batch_controls:
        for dose in control_dict[control].keys():
            z[0] += control_dict[control][dose]
            z[1] += get_well_id.get_well_id(reps[0], control, float(dose))
            z[2] += get_well_id.get_well_id(reps[1], control, float(dose))
            z[3] += get_well_id.get_well_id(reps[2], control, float(dose))
    return z


control_dict = {}
control_dict['DMSO'] = {'300': ['A11', 'B11', 'C11', 'D11', 'E11', 'F11', 'G11', 'H11', 
                                'A12', 'B12', 'C12', 'D12', 'E12', 'F12', 'G12', 'H12',
                                'A10', 'B10', 'C10', 'D10']}
control_dict['pc_Paclitaxel'] = {'1': ['A11', 'B11', 'C11', 'D11']}
control_dict['pc_GSK2126458'] = {'0.1': ['A12', 'B12', 'C12', 'D12']}
control_dict['pc_Actinomycin_D'] = {'1': [chr(65+r) + '11'
                                          for r in range(8)]}
control_dict['pc_Staurosporine'] = {'1': [chr(65+r) + '12'
                                          for r in range(8)]}


batch1_drugs = drugs[:8]
batch2_drugs = drugs[8:16]
batch3_drugs = drugs[16:24]
batch4_drugs = drugs[24:]
batch1_controls = ['pc_Paclitaxel', 'pc_GSK2126458']
batch2_controls = ['DMSO']
batch3_controls = ['pc_Actinomycin_D', 'pc_Staurosporine']



def get_map_report(drugs, Designs, drug_treatment,
                   control_dict, file_prefix, batch_controls=[]):
    
    z = get_mapping(drugs, Designs, drug_treatment)
    if batch_controls:
        z = get_control_mapping(batch_controls, control_dict, z)


    columns = ['Source Plate', 'Source Copy', 'Source Well', 'Source Plate type',
               'Destination Well', 'Destination Plate type', 
               'Person Visiting', 'Screen Number', 'Volume']
    sp = ['']*len(z[1])
    sc = ['']*len(z[1])
    person = ['LD']*len(z[1])
    spt = ['Nunc 96 UB PP']*len(z[1])
    dpt = ['']*len(z[1])
    screen_number = ['']*len(z[1])
    volume = [46]*len(z[1])
    df1 = pd.DataFrame(zip(sp, sc, z[0], spt, z[1], dpt,
                           person, screen_number, volume), columns=columns)
    df1.to_csv('OUTPUT/%s_DUBI1_A.csv' % file_prefix, index=False)
    df2 = pd.DataFrame(zip(sp, sc, z[0], spt, z[2], dpt,
                           person, screen_number, volume), columns=columns)
    df2.to_csv('OUTPUT/%s_DUBI1_B.csv' % file_prefix, index=False)
    df3 = pd.DataFrame(zip(sp, sc, z[0], spt, z[3], dpt,
                           person, screen_number, volume), columns=columns)
    df3.to_csv('OUTPUT/%s_DUBI1_C.csv' % file_prefix, index=False)


get_map_report(batch1_drugs, Designs, drug_treatment,
               control_dict, 'SP1', batch1_controls)
get_map_report(batch2_drugs, Designs, drug_treatment,
               control_dict, 'SP2', batch2_controls)
get_map_report(batch3_drugs, Designs, drug_treatment,
               control_dict, 'SP3', batch3_controls)
get_map_report(batch4_drugs, Designs, drug_treatment,
               control_dict, 'SP4')

