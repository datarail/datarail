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


def generate_GRinput(df, counts_column='cell_count',
                     time0_plate=True, mean=True):
    dft = df[df['role'] == 'treatment'].copy()
    df_counts = dft.copy()
    cell_lines = dft['cell_line'].unique()
    if time0_plate:
        time0_dc = {}
        for line in cell_lines:
            dfc = df[df['cell_line'] == line].copy()
            time0_ctrl = dfc[dfc['timepoint'] == 'time0_ctrl'][
                counts_column].mean()
            time0_dc[line] = time0_ctrl
        df_counts['cell_count__time0'] = [time0_dc[line] for line in
                                          df_counts['cell_line'].tolist()]
        counts_col = ['cell_count', 'cell_count__ctrl', 'cell_count__time0']
    else:
        counts_col = ['cell_count', 'cell_count__ctrl']
    dc2 = df[df.role == 'negative_control'].groupby(
        ['cell_line', 'timepoint'])['cell_count'].mean()
    df_counts['cell_count__ctrl'] = [dc2.loc[c, t] for c, t in
                                     zip(df_counts['cell_line'],
                                         df_counts['timepoint'])]
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

    if mean:
        df_counts_mean = df_counts.groupby(['cell_line', 'agent',
                                            'concentration', 'timepoint'])[
                                                counts_col].apply(np.mean)
        df_counts_mean = df_counts_mean.reset_index()
        return df_counts_mean
    else:
        return df_counts


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
                      dead_cell_columns,
                      total_cell_columns):
    """ Appends column for fraction of dead cells
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
