import pandas as pd
import numpy as np
from collections import OrderedDict
from itertools import groupby
import edge_barcode
import warnings


def make_assay(csv_file, plate_dims, barcode_prefix,
               encode_plate=False, num_replicates=1):
    df = pd.read_csv(csv_file)
    drugs = df.Compound_Name[df.Role == 'treatment'].tolist()
    positive_controls = df.Compound_Name[df.Role ==
                                         'positive_control'].tolist()
    negative_controls = df.Compound_Name[df.Role ==
                                         'negative_control'].tolist()
    barcode_treatments = df.Compound_Name[df.Role == 'barcode'].tolist()

    drug_treatments = OrderedDict()
    for drug in drugs:
        max_dose = df['Highest_Dose'].ix[
            df['Compound_Name'] == drug].values[0]
        num_doses = df['num_wells'].ix[
            df['Compound_Name'] == drug].values[0]
        max_dose_value, _ = split_text(max_dose)
        drug_treatments[drug] = float(max_dose_value) * 1e-4 * np.logspace(
            0, 4, num_doses)

    num_wells = plate_dims[0] * plate_dims[1]
    num_treatments = sum(len(v) for v in drug_treatments.itervalues())
    num_edge_wells = get_boundary_cell_count(plate_dims)
    inner_wells_available = num_wells - num_treatments - num_edge_wells

    nc_treatments = OrderedDict()
    for nc in negative_controls:
        max_dose = df['Highest_Dose'].ix[
            df['Compound_Name'] == nc].values[0]
        num_wells = df['num_wells'].ix[
            df['Compound_Name'] == nc].values[0]
        max_dose_value, _ = split_text(max_dose)
        nc_treatments[nc] = [max_dose_value] * num_wells
    num_nc_treatments = sum(len(v) for v in nc_treatments.itervalues())
    total_control_wells = num_nc_treatments

    if num_nc_treatments < 8:
        print ""
        warnings.warn(
            'Insufficent number of wells alloted for negative controls')
        print "Atleast 8 wells have to be assigned for negative controls,"\
            " recommended number is 12, user has currently alloted %d wells"\
            " for negative_controls" % num_nc_treatments

    pc_treatments = OrderedDict()
    try:
        for pc in positive_controls:
            max_dose = df['Highest_Dose'].ix[
                df['Compound_Name'] == pc].values[0]
            num_wells = df['num_wells'].ix[
                df['Compound_Name'] == pc].values[0]
            max_dose_value, _ = split_text(max_dose)
            pc_treatments[pc] = [max_dose_value] * num_wells
        num_pc_treatments = sum(len(v) for v in pc_treatments.itervalues())
        total_control_wells += num_pc_treatments
    except NameError:
        pass

    if total_control_wells > inner_wells_available:
        print ""
        warnings.warn(
            'Number of wells alloted for controls exceeds available wells')
        print "%d wells are available for controls, user has alloted %d wells"\
            " for negative controls and %d for positive controls" % (
                inner_wells_available, num_nc_treatments, num_pc_treatments)

    if encode_plate:
        bc_treatments = OrderedDict()
        try:
            for bc in barcode_treatments:
                max_dose = df['Highest_Dose'].ix[
                    df['Compound_Name'] == bc].values[0]
                num_wells = df['num_wells'].ix[
                    df['Compound_Name'] == bc].values[0]
                max_dose_value, _ = split_text(max_dose)
                bc_treatments[bc] = [max_dose_value] * num_wells
            num_bc_treatments = sum(len(v) for v in bc_treatments.itervalues())

            barcodes = [barcode_predix + chr(65+i)
                        for i in range(num_replicates)]
            num_bc_wells = [edge_barcode.encode_barcode(bc) for bc in barcodes]
            max_bc_wells = n.max(num_bc_wells)
            if num_bc_treatments < max_bc_wells:
                print ""
                warnings.warn(
                    'Insufficent number of wells alloted for encodng barcode')
                print "barcode requires %d wells, user has alloted %d wells"\
                    % (max_bc_wells, num_bc_treatments)
        except NameError:
            print ""
            warnings.warn('treatments for barcode not specified')

    return drug_treatments, nc_treatments, pc_treatments, bc_treatments


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


def split_text(s):    
    for k, g in groupby(s, str.isalpha):
        yield ''.join(list(g))


def get_boundary_cell_count(plate_dims):
    """ get number of wells in outer or inner edges

    Parameter
    ---------
    plate_dims: array
         dimensions of plate

    Returns
    -------
    boundary_cell_count: int
           number of wells in the edges
    """
    boundary_cell_count = 2 * (plate_dims[0] + plate_dims[1] - 2)
    return boundary_cell_count
