import pandas as pd
import collections
import numpy as np


def merge_plate_level_metadata(dfo, dfp, common_identifier='barcode'):
    dfo['well'] = ["%s%s" % (w[0], w[1:].zfill(2))
                   for w in dfo.well.tolist()]
    dfp.index = dfp[common_identifier].tolist()
    dfo.index = dfo[common_identifier].tolist()
    dfp_lt = dfp.loc[dfo.index.tolist()]
    dfc = pd.concat([dfo, dfp_lt], axis=1)
    dfc = dfc.loc[:, ~dfc.columns.duplicated()]
    return dfc


# def merge_well_level_metadata(dfc, dfmeta):
#     dfm = dfmeta.copy()
#     dfm.index = ["%s_%s" % (r, w) for r, w in
#                  zip(dfm['plate'], dfm['well'])]
#     del dfm['plate']
#     df_list = []
#     for barcode in dfc.barcode.unique():
#         dfcb = dfc[dfc.barcode == barcode].copy()
#         dfcb.index = ["%s_%s" % (r, w) for r, w in
#                       zip(dfcb['plate'], dfcb.well)]
#         try:
#             dfcb = remove_duplicate_wells(dfcb)
#             print('Removing duplicate wells in %s' % barcode)
#         except KeyError:
#             pass
#         df_list.append(pd.concat([dfcb, dfm], axis=1))
#     dfc2 = pd.concat(df_list)
#     dfc3 = dfc2.dropna(subset=['barcode'])
#     return dfc3


def merge_well_level_metadata(dfc, dfmeta):
    dfc['well'] = ["%s%s" % (s[0], s[1:].zfill(2))
                   for s in dfc.well.tolist()]
    if 'Result' in dfc.columns.tolist():
        dfc[['barcode', 'timestamp']] = dfc.Result.str.split(
            ' > ', expand=True)
    dfm = dfmeta.copy()
    dfm.index = ["%s_%s" % (r, w) for r, w in
                 zip(dfm['barcode'], dfm['well'])]
    dfc.index = ["%s_%s" % (r, w) for r, w in
                 zip(dfc['barcode'], dfc['well'])]
    dfc2 = pd.concat([dfc, dfm], axis=1, sort=False)
    dfc3 = dfc2.loc[:, ~dfc2.columns.duplicated()].copy()
    # dfc3 = dfc2.dropna(subset=['Result'])
    return dfc3


def remove_duplicate_wells(dfo):
    # find duplicate_wells
    dup_wells = [item for item, count in
                 collections.Counter(dfo.well.tolist()).items()
                 if count > 1]
    dfod = dfo[dfo.well.isin(dup_wells)]
    dup_urls = dfod.sort_values(['Result']).groupby(['well']).apply(
        lambda group: group.iloc[0, 1:])['URL'].tolist()
    dfo = dfo[~dfo.URL.isin(dup_urls)].copy()
    return dfo


def generate_GRinput(df_input, counts_column='cell_count',
                     time0_plate=True, 
                     evaluate_dead=True,
                     dead_cell_columns=['corpse_count', 'cell_count__dead']):
    """Takes cell/dead count values across all plates to compute mean counts 
    for time0_ctrl and untreated controls. The columns are then labeled as per 
    expected input for computing gr values and metrics

    Parameters
    ----------
    df : pandas dataframe
      input dataframe consisting of well level live and dead cell counts.
    live_counts_column : str
      column name corresponding to number of live cells.
    time0_plate : boolean
      Specify as True if time_0 plate is included in the experiment.Default is True.
    evaluate_dead : boolean
      computes mean dead cell count for time0_ctrl and untreated controls. Default is True.
    dead_cell_columns : list of str
        list of columns corresponding to dead cell counts.

    Returns
    -------
    df_counts : pandas dataframe
        live/dead cell counts per condition with corresponding time0_ctrl and untreated
        control counts. Serves as input for computing GR values and metrics.
    """
    df = df_input.copy()
    df['timepoint'] = df['timepoint'].replace([0], 'time0_ctrl')
    df.loc[df.timepoint == 'time0_ctrl','role'] = 'negative_control'
    if evaluate_dead:
        df['dead_count'] = df[dead_cell_columns].sum(axis=1)
        total_cell_columns = dead_cell_columns + [counts_column]
        df = get_fraction_dead(df, dead_cell_columns, total_cell_columns).copy()
    dft = df[df['role'] == 'treatment'].copy()
    df_counts = dft.copy()
    cell_lines = dft['cell_line'].unique()
    if time0_plate:
        time0_dict_live = {}
        time0_dict_dead = {}
        for line in cell_lines:
            dfc = df[df['cell_line'] == line].copy()
            time0_ctrl_live = dfc[dfc['timepoint'] == 'time0_ctrl'][
                counts_column].mean()
            time0_dict_live[line] = time0_ctrl_live
            if evaluate_dead:
                time0_ctrl_dead = dfc[dfc['timepoint'] == 'time0_ctrl'][
                    'dead_count'].mean()
                time0_dict_dead[line] = time0_ctrl_dead
        df_counts['cell_count__time0'] = [time0_dict_live[line] for line in
                                          df_counts['cell_line'].tolist()]
        if evaluate_dead:
            df_counts['dead_count__time0'] = [time0_dict_dead[line] for line in
                                          df_counts['cell_line'].tolist()]
        #counts_col = ['cell_count', 'cell_count__ctrl', 'cell_count__time0']
    #else:
    #    counts_col = ['cell_count', 'cell_count__ctrl']
    dc2 = df[df.role == 'negative_control'].groupby(
        ['barcode', 'cell_line', 'timepoint'])['cell_count'].mean()
    df_counts['cell_count__ctrl'] = [dc2.loc[b, c, t] for b, c, t in
                                     zip(df_counts['barcode'],
                                         df_counts['cell_line'],
                                         df_counts['timepoint'])]
    if evaluate_dead:
        dcd = df[df.role == 'negative_control'].groupby(
            ['barcode', 'cell_line', 'timepoint'])['dead_count'].mean()
        df_counts['dead_count__ctrl'] = [dcd.loc[b, c, t] for b, c, t in
                                         zip(df_counts['barcode'],
                                             df_counts['cell_line'],
                                             df_counts['timepoint'])]
        dcf = df[df.role == 'negative_control'].groupby(
            ['barcode', 'cell_line', 'timepoint'])['fraction_dead'].mean()
        df_counts['fraction_dead__ctrl'] = [dcf.loc[b, c, t] for b, c, t in
                                            zip(df_counts['barcode'],
                                                df_counts['cell_line'],
                                                df_counts['timepoint'])]
        df_counts['increase_fraction_dead'] = df_counts['fraction_dead'] -\
            df_counts['fraction_dead__ctrl']
    df_counts = df_counts.loc[:, ~df_counts.columns.duplicated()]

    agent_cols = [a for a in df_counts.columns.tolist()
                  if a.startswith('agent')]
    concentration_cols = [c for c in df_counts.columns.tolist()
                          if c.startswith('concentration')]
    df_counts[agent_cols] = df_counts[agent_cols].where(
        pd.notnull(df_counts), '')
    df_counts[concentration_cols] = df_counts[concentration_cols].where(
        pd.notnull(df_counts), 0)
    df_counts['replicate'] = df_counts.groupby(['cell_line'] + agent_cols +
                                               concentration_cols)[
        'cell_count'].transform(lambda x: pd.factorize(x)[0]+1)
    df_counts['replicate'] = df_counts['replicate'].replace([0], np.nan)
    if len(agent_cols) > 1:
        df_counts['agent'] = df_counts[agent_cols].apply(lambda x: '; '.join(
            x[x.notnull()]), axis=1)
        df_counts['concentration'] = df_counts[concentration_cols].astype(str).apply(
            lambda x: '; '.join(x[x.notnull()]), axis=1)
    df_counts = df_counts.sort_values(by=['cell_line', 'agent', 'concentration'])
    df_counts['cell_count'] = df_counts['cell_count'].fillna(0)
    df_counts = df_counts[df_counts.cell_count != 0].copy() # TEMP HASH to be checked
    return df_counts    


def get_counts_mean(df_counts):
    scol = ['cell_count', 'cell_count__ctrl', 'cell_count__time0',
            'dead_count', 'dead_count__ctrl', 'dead_count__time0',
            'fraction_dead', 'fraction_dead__ctrl', 'increase_fraction_dead']
    counts_col = [s for s in scol if s in df_counts.columns.tolist()]
    df_counts_mean = df_counts.groupby(['cell_line', 'agent',
                                        'concentration', 'timepoint'])[
                                                counts_col].apply(np.mean)
    df_counts_mean = df_counts_mean.reset_index()
    return df_counts_mean
    


def make_long_table(output_file, barcode, plate_dims=[16, 24]):
    """Use with readouts that are in plate layout format
    eg:- CTG
    """
    df = pd.read_excel(output_file)
    df = df.iloc[:plate_dims[0], :plate_dims[1]]
    columns = df.columns.tolist()
    df['Row'] = df.index.tolist()
    dfp = pd.melt(df, id_vars=['Row'], value_vars=columns)
    dfp['well'] = ["%s%02d" % (r, c) for r, c
                   in zip(dfp.Row.tolist(), dfp.variable.tolist())]
    dfp = dfp.rename(columns={'value': 'cell_count'})
    dfp = dfp[['well', 'cell_count']]
    dfp['barcode'] = [barcode]*len(dfp)
    return dfp


def get_fraction_dead(dfo,
                      dead_cell_columns=['corpse_count', 'cell_count__dead'], 
                      total_cell_columns=['cell_count', 'cell_count__dead',
                                          'corpse_count']):
    """Appends column for fraction of dead cells

    Parameters
    ----------
    dfo : pandas dataframe
    dead_cell_columns : list of str
        list of column names corresponding to dead cells
    total_cell_columns : list of str
        list of column names corresponding to total cell counts

    Returns
    -------
    dfo2 : pandas dataframe    
    """
    dfo2 = dfo.copy()
    dfo2['fraction_dead'] = dfo2[dead_cell_columns].sum(axis=1).div(
        dfo2[total_cell_columns].sum(axis=1))
    return dfo2


def get_time0_metadata(dfm, time0_barcode, ref_barcode):
    """Construct time0 metadata table from a reference treated metadata table
    
    Parameters
    ----------
    dfm : pandas.DataFrame
       Metadata from datarail or D300 software.
    time0_barcode : str
       Barcode of time0 plate for which metadata table needs to be constructed.
    ref_barcode : str
       Barocode or plate identifier of reference plate.

    Return
    ------
    dt : pandas.DataFrame
      time0 metadata table 
    """
    dt = dfm[dfm.barcode == ref_barcode].copy()
    dt['timepoint'] = 'time0_ctrl'
    dt['barcode'] = time0_barcode
    treated_wells = dt[dt.cell_line.notnull()].well.unique()
    dt['role'] = np.nan
    dt['agent'] = np.nan
    dt['concentration'] = 0
    dt.loc[da.well.isin(treated_wells), 'role'] = 'negative_control'
    dt.loc[da.well.isin(treated_wells), 'agent'] = np.nan
    return dt
    
