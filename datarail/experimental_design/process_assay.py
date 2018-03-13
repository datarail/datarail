import pandas as pd
import numpy as np
from collections import OrderedDict
from itertools import groupby
from itertools import product
import datarail.experimental_design.edge_fingerprint as edge_fingerprint
import warnings


def get_boundary_cell_count(plate_dims, exclude_outer=1):
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
    if exclude_outer == 2:
        boundary_cell_count += 2 * (plate_dims[0]-2 + plate_dims[1]-2 - 2)
    return boundary_cell_count


def set_dosing(max_dose, num_doses, num_replicates=1):
    """ returns list of doses (micromolar) at half log intervals 
    Parameters
    ----------
    max_dose: int
        highest dose in the dose range
    num_doses: int
        number of doses in the dose range. Maximum is set to 9. 
        If num doses < 9, then doses will be removed starting from lowest moving to higher doses.
    num_replicates: int
        number of times the dose range is replicated on the plate
    
    Returns:
    -------
    dose_range: list of floats
    """
    dose_range = max_dose * 1e-4 * np.logspace(0, 4, 9)
    dose_range = [round(s, 4) for s in dose_range]
    dose_range = sorted(list(set(dose_range)))[::-1]
    dose_range = dose_range[:num_doses]
    if num_replicates > 1:
        dr = list(set(dose_range))
        for i in range(1, num_replicates):
            dose_range.extend(dr)
    return dose_range


def set_combo_dosing(max_doses, num_doses, eq=False, num_replicates=1):
    """ returns combination of doses  for 2 or more agents 
    Parameters
    ----------
    max_doses: list of int
        highest doses in the dose range for each agent
    num_doses: list of int
        number of doses in the dose range for each agent. Maximum is set to 9. 
        If num doses < 9, then doses will be removed starting from lowest moving to higher doses.
    num_replicates: list of int
        list of number of times the dose range is replicated on the plate for each agent
    
    Returns:
    -------
    dose_range: list of tuples
        each tuple corresponds to one combination of doses for 2 or more agents
    """
    dose_lists = []
    for md, nd in zip(max_doses, num_doses):
        dose_range = md * 1e-4 * np.logspace(0, 4, 9)
        dose_range = sorted(list(set(dose_range)))[::-1]
        dose_range = [round(s, 4) for s in dose_range]
        dose_range = dose_range[:nd]
        dose_lists.append(dose_range)        
    combo_doses = list(product(*dose_lists))
    if num_replicates > 1:
        cd = list(set(combo_doses))
        for i in range(1, num_replicates):
            combo_doses.extend(cd)
    if eq:
        combo_doses = [c for c in combo_doses if len(set(c)) == 1]
    return combo_doses    


def exclude_treatment(df, drug, doses):
    df2 = df.copy()
    df2 = df2[~((df2.agent == drug) & (df2.concentration.isin(doses)))]
    return df2


def construct_well_level_df(input_file, plate_dims=[16, 24],
                            exclude_outer=1):
    """ Generates long table of doses and drugs mapped to wells.

    Input file should be broad description of the planned experimental design
    and should contain the following columns
    - 'agent' column listing the names of drugs
       including positve and negative controls
    - 'max_dose__um' column listing the highest dose for each agent
    - 'num_doses' column listing the number of doses for each agent
    - 'role' column listing the intended role for each agent
       'treatment', positive_control', 'negative_control', or 'fingerprint'
    - 'num_replicates' column listing number of times a drug's
       dosing scheme is replicated on the same plate
    - 'exclude_doses' column (OPTIOINAL) listing doses to be excluded for a given drug

    Parameters:
    ----------
    input_file: str
            csv or tsv file name of the input file
    plate_dims: list
            dimensions of the physical plate
    exlude_outer: int
            number of outer wells to exlude; defaut set to 1
    Returns
    -------
    df_well: pandas dataframe
           dataframe mapping wells to drug and dose response
    """
    df_spec = pd.read_csv(input_file)
    df_spec = df_spec.fillna('')
    if 'exclude_doses' not in df_spec.columns.tolist():
        df_spec['exclude_doses'] = ['']*len(df_spec)
    df_spec['exclude_doses'] = df_spec['exclude_doses'].fillna('')
    drugs, doses, role, identifier = [], [], [], []
    df_tr = df_spec[df_spec.role == 'treatment'].copy()
    for drug in df_tr.agent.tolist():
        max_dose = df_tr[df_tr.agent == drug]['max_dose__um'].values[0]
        max_dose = str(max_dose).split(',')
        max_dose = [float(mx) for mx in max_dose]
        num_doses = df_tr[df_tr.agent == drug]['num_doses'].values[0]
        num_doses = str(num_doses).split(',')
        num_doses = [int(nd) for nd in num_doses]
        num_replicates = df_tr[df_tr.agent == drug][
            'num_replicates'].values[0]
        print(num_replicates)
        exclude_doses = df_tr[df_tr.agent == drug][
            'exclude_doses'].values[0]
        if len(max_dose) == 1:
            max_dose = max_dose[0]
            num_doses = num_doses[0]
            dose_range = set_dosing(max_dose, num_doses, num_replicates)
            #if exclude_doses == '':
            #    dose_range = dose_range[:num_doses]
            #else:
            #    exclude_doses = [float(s) for s in exclude_doses.split(',')]
            #    dose_range = [d for d in dose_range if d not in exclude_doses]
        else:
            eqm = df_tr[df_tr.agent == drug]['equivalent'].values[0]
            dose_range = set_combo_dosing(max_dose, num_doses,
                                          eq=bool(eqm),
                                          num_replicates=num_replicates)
        doses += dose_range
        drugs += [drug] * len(dose_range)
        role += ['treatment'] * len(dose_range)
        identifier += ["%s_%d" % (drug, num) for num in range(len(dose_range))]
    if 'positive_control' in df_spec.role.unique():
        dfp = df_spec[df_spec.role == 'positive_control'].copy()
        for drug in dfp.agent.tolist():
            max_dose = dfp[dfp.agent == drug]['max_dose__um'].values[0]
            num_replicates = int(dfp[dfp.agent == drug][
                'num_replicates'].values[0])
            doses += [max_dose] * num_replicates
            drugs += [drug] * num_replicates
            role += ['positive_control'] * num_replicates
            identifier += ["%s_%d" % (drug, num)
                           for num in range(num_replicates)]
    else:
        warnings.warn(
            'Experimental design does not have positive_controls')
    num_outer_wells = get_boundary_cell_count(plate_dims, exclude_outer)
    num_available_wells = (plate_dims[0] * plate_dims[1]) - num_outer_wells
    num_treatment_wells = len(doses)
    if num_available_wells < num_treatment_wells:
        warnings.warn('Number of treatment wells required (%d)'
                      'exceed available wells (%d)' % (
                         num_treatment_wells, num_available_wells))
    df_well = pd.DataFrame(list(zip(drugs, doses, role, identifier)),
                           columns=['agent', 'concentration',
                                    'role', 'identifier'])
    return df_well


def add_negative_control(df, control_name='DMSO',
                         plate_dims=[16, 24], exclude_outer=1):
    """ Assigns negative control agent to untreated wells
    Parameter
    ---------
    df: pandas dataframe
        well level metadata with specification of agents and concentration
    control_name: str
        name of control agent; default is DMSO
    plate_dims: list of int
        dimension of physical plate
    exclude_outer: int
       number of outer well layers to be excluded
    
    Returns:
    -------
    df_well: pandas dataframe
        well level metadata with specification of both agents and negative control
    """
    num_treatment_wells = len(df)
    num_outer_wells = get_boundary_cell_count(plate_dims, exclude_outer)
    num_available_wells = (plate_dims[0] * plate_dims[1]) - num_outer_wells
    num_nc_wells = num_available_wells - num_treatment_wells
    if num_nc_wells < 8:
        print("")
        warnings.warn(
            'Insufficent number of wells alloted for negative controls')
        print("Atleast 8 wells have to be assigned for negative controls,"
              " recommended number is 12, user has currently alloted %d wells"
              " for negative_controls" % num_nc_wells)
    role = df.role.tolist()
    doses = df.concentration.tolist()
    drugs = df.agent.tolist()
    identifiers = df.identifier.tolist()
    role += ['negative_control'] * num_nc_wells
    doses += [np.nan] * num_nc_wells
    drugs += [control_name] * num_nc_wells
    identifiers += [control_name] * num_nc_wells
    df_well = pd.DataFrame(list(zip(drugs, doses, role, identifiers)),
                           columns=['agent', 'concentration', 'role',
                                    'identifier'])
    return df_well


def assign_fingerprint_wells(fingerprint, treatment, dose):
    """ Returns a set of wells along the edge that serve as barcode for
    the plate based on the fingerprint word
    Parameters
    ----------
    fingerprint: str
    treatment: str
        the drug used to treat fingerprint wells
    dose: float
        the dose of drug treatment
    Returns:
    -------
    df: pandas dataframe
       table mapping wells encoding the fingerprint to treatment
    """
    fingerprint_wells = edge_fingerprint.encode_fingerprint(fingerprint)
    treatment_list = [treatment] * len(fingerprint_wells)
    dose_list = [dose] * len(fingerprint_wells)
    role = ['fingerprint'] * len(fingerprint_wells)
    identifier = ["%s_%d" % (treatment, d)
                  for d in range(len(fingerprint_wells))]
    df = pd.DataFrame(list(zip(fingerprint_wells, treatment_list,
                               dose_list, role, identifier)),
                      columns=['well', 'agent', 'concentration',
                               'role', 'identifier'])
    return df


def define_treatment_wells(exclude_outer=1):
    """ defines set of inner wells to be used for treatments
    Parameter:
    ---------
    exclude_outer: int
       defines outer well columns and rows to to exclude
    Returns:
    -------
    tr_wells, list(set(exclude_wells)): tuple of lists
       lists of treatment wells and outer wells
    """
    cols = ["%02d" % s for s in range(1, 25)]
    rows = [chr(65+n) for n in range(16)]
    wells = []
    for row in rows:
        for col in cols:
            wells.append("%s%s" % (row, col))
    exclude_pttrn = ['A', 'P', '01', '24']
    if exclude_outer == 2:
        exclude_pttrn = ['A', 'B', 'O', 'P', '01', '02', '23', '24']
    exclude_wells = []
    for well in wells:
        for ep in exclude_pttrn:
            if well.find(ep) != -1:
                exclude_wells.append(well)
    tr_wells = [s for s in wells if s not in exclude_wells]
    return tr_wells, list(set(exclude_wells))


def randomize_wells(df, plate_names=['BCA2_A'],
                    fingerprint_drug=None, fingerprint_dose=1, exclude_outer=1):
    """ Returns dataframe with randomized wells for all plate replicates
    Parameters:
    -----------
    df: pandas dataframe
        ordered wells mapped to drug and dose for each master plate
    fingerprints: list of str
        list of plate identifiers (barcode) for all plate replicates
    fingerprint_drug: str
        drug used for treating fingerprint wells
    fingerprint_dose: int
       dose of drug used for fingerprint wells treatment
    Returns:
    --------
    dfr: pandas dataframe
       drug and dose mapped to randomized wells
    """
    tr_wells, _ = define_treatment_wells(exclude_outer=exclude_outer)
    df['well'] = tr_wells
    ordered_wells = df.well.tolist()
    cols = ["%02d" % s for s in range(1, 25)]
    rows = [chr(65+n) for n in range(16)]
    wells = []
    for row in rows:
        for col in cols:
            wells.append("%s%s" % (row, col))
    df_list = []
    for rep_num, plate in enumerate(plate_names):
        np.random.seed(rep_num+1)
        randomized_wells = np.random.choice(ordered_wells,
                                            size=len(ordered_wells),
                                            replace=False)
        df['well'] = randomized_wells
        df.index = df['well']
        if fingerprint_drug:
            df_fp = assign_fingerprint_wells(plate,
                                             fingerprint_drug,
                                             fingerprint_dose)
            df_fp.index = df_fp['well']
            dfc = pd.concat([df_fp, df])
        else:
            dfc = df.copy()
        dfc['plate'] = [plate] * len(dfc)
        remainder_wells = [w for w in wells if w not in dfc.well.tolist()]
        dfo = pd.DataFrame(list(zip(remainder_wells, [plate] * len(
            remainder_wells))), columns=['well', 'plate'])
        dfc2 = pd.concat([dfo, dfc])
        dfc2 = dfc2.sort_values(['well'])
        df_list.append(dfc2)
    dfr = pd.concat(df_list)
    dfr[['agent', 'role']] = dfr[['agent', 'role']].where(pd.notnull(dfr), '')
    dfr['concentration'] = dfr['concentration'].fillna(0)
    dfr.index = range(len(dfr))
    max_agents = np.max([len(a.split(',')) for a in dfr.agent.tolist()])
    if max_agents > 1:
        agent_columns = ['agent%d' % ma for ma in range(1, max_agents+1)]
        dfr['agent'] = dfr['agent'].str.replace(' ', '')
        dfr[agent_columns] = dfr.agent.str.split(',', expand=True)
        del dfr['agent']
        concentration_columns = ['concentration%d' % ma for ma in range(1, max_agents+1)]
        dfr[concentration_columns] = dfr['concentration'].apply(pd.Series)
        del dfr['concentration']
    return dfr


# def split_text(s):
#     for k, g in groupby(s, str.isalpha):
#         yield ''.join(list(g))
