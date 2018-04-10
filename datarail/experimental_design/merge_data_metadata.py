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
    dfm = dfmeta.copy()
    dfm.index = ["%s_%s" % (r, w) for r, w in
                 zip(dfm['barcode'], dfm['well'])]
    dfc.index = ["%s_%s" % (r, w) for r, w in
                 zip(dfc['barcode'], dfc['well'])]
    dfc2 = pd.concat([dfc, dfm], axis=1)
    dfc3 = dfc2.dropna(subset=['Result'])
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


def generate_GRinput(df, counts_column='cell_count'):
    dft = df[df['role'] == 'treatment'].copy()
    df_counts = dft.copy()
    cell_lines = dft['cell_line'].unique()
    dc1 = {}
    for line in cell_lines:
        dfc = df[df['cell_line'] == line].copy()
        time0_ctrl = dfc[dfc['timepoint'] == 'time0_ctrl'][
            counts_column].mean()
        dc1[line] = time0_ctrl
    dc2 = df[df.role == 'negative_control'].groupby(
        ['cell_line', 'timepoint'])['cell_count'].mean()
    df_counts['cell_count__ctrl'] = [dc2.loc[c, t] for c, t in
                                     zip(df_counts['cell_line'],
                                         df_counts['timepoint'])]
    df_counts['cell_count__time0'] = [dc1[line] for line in
                                      df_counts['cell_line'].tolist()]
    df_counts = df_counts.loc[:, ~df_counts.columns.duplicated()]
    df_counts_mean = df_counts.groupby(['cell_line', 'agent',
                                        'concentration', 'timepoint'])[
        ['cell_count__ctrl', 'cell_count__time0', 'cell_count']].apply(np.mean)
    df_counts_mean = df_counts_mean.reset_index()
    return df_counts_mean
