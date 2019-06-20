import pandas as pd
import numpy as np


def _combine_combos(df):
    dfc = df.groupby(['barcode', 'well'], as_index=True).apply(
          lambda x: x.apply(_combine_duplicates, axis=0))
    return dfc


def _combine_duplicates(x):
    return ','.join(x.astype(str))


def export2pd(filename):
    """Note: Open .xml file in excel and save as .xlsx file

    Parameters
    ----------
    filename : str
        path to .xlsx file

    Returns
    -------
    dfm : pandas datafame
       metadata file complying with datarail metadata conventions
    """
    dfi = pd.read_excel(filename, sheet_name='Tabular detail')
    dfi = dfi.rename(columns={'Dispensed concentration': 'Dispensed conc.'})
    dfi = dfi[['Plate', 'Dispensed\nwell', 'Fluid name', 'Dispensed conc.']]
    dfi.columns = ['barcode', 'well', 'agent', 'concentration']
    dfi['concentration'] = dfi['concentration'].fillna(0)
    dfm = _combine_combos(dfi)
    del dfm['barcode']
    del dfm['well']
    dfm = dfm.reset_index()
    dfm['agent'] = [s.replace(',DMSO normalization', '')
                    for s in dfm.agent.tolist()]
    dfm['agent'] = [s.replace(',a+Tw normalization', '')
                    for s in dfm.agent.tolist()]
    dfm['concentration'] = [','.join(c.split(',')[:len(a.split(','))])
                            for a, c  in zip(dfm.agent, dfm.concentration)]
    
    dfm['agent'] = [s.replace('DMSO normalization', 'DMSO')
                    for s in dfm.agent.tolist()]
    dfm['role'] = ['negative_control' if s == 'DMSO' else 'treatment'
                   for s in dfm.agent.tolist()]
    max_agents = np.max([len(a.split(',')) for a in dfm.agent.tolist()])
    if max_agents > 1:
        agent_columns = ['agent%d' % ma for ma in range(1, max_agents+1)]
        dfm['agent'] = dfm['agent'].str.replace(' ', '')
        dfm[agent_columns] = dfm.agent.str.split(',', expand=True)
        dfm[agent_columns] = dfm[agent_columns].where(pd.notnull(dfm), '')
        del dfm['agent']
        concentration_columns = ['concentration%d' % ma
                                 for ma in range(1, max_agents+1)]
        dfm[concentration_columns] = dfm.concentration.str.split(',',
                                                                 expand=True)
        dfm[concentration_columns] = dfm[concentration_columns].fillna(0)
        del dfm['concentration']
        for col in concentration_columns:
            dfm[col] = [round(float(s), 4) for s in dfm[col].tolist()]
    else:
        dfm['concentration'] = [round(float(s), 4)
                                for s in dfm.concentration.tolist()]
    dfp = pd.read_excel(filename, sheet_name='Plates')
    plate_map = {pn: pb for pn, pb in zip(dfp['Plate'],
                                          dfp['Plate name'])}
    dfm['barcode'] = dfm['barcode'].map(plate_map)
    dfm['agent'] = dfm[agent_columns].apply(
        lambda x: ','.join(x.dropna().astype(str)), axis=1)
    dfm[concentration_columns] = dfm[concentration_columns].replace([0], np.nan)
    dfm['concentration'] = dfm[concentration_columns].apply(
        lambda x: ','.join(x.dropna().astype(str)), axis=1)
    return dfm
