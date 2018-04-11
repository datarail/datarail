import pandas as pd
import numpy as np


def combine_combos(df):
    dfc = df.groupby(['barcode', 'well'], as_index=True).apply(
          lambda x: x.apply(combine_duplicates, axis=0))
    return dfc


def combine_duplicates(x):
    return ';'.join(x.astype(str))


def export2pd(filename):
    """
    Note: Open .xml file in excel and save as .xlsx file
    Parameter
    ---------
    filename: str
        path to .xlsx file
    Returns
    -------
    dfm: pandas datafame
       metadata file complying with datarail metadata conventions
    """
    dfi = pd.read_excel(filename, sheet_name='Tabular detail')
    dfi = dfi[['Plate', 'Dispensed\nwell', 'Fluid name', 'Dispensed conc.']]
    dfi.columns = ['barcode', 'well', 'agent', 'concentration']
    dfi['concentration'] = dfi['concentration'].fillna(0)
    dfm = combine_combos(dfi)
    del dfm['barcode']
    del dfm['well']
    dfm = dfm.reset_index()
    dfm['agent'] = [s.replace(';DMSO normalization', '')
                    for s in dfm.agent.tolist()]
    dfm['concentration'] = [s.replace(';0.0', '') for s in
                            dfm.concentration.tolist()]
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
        dfm[concentration_columns] = dfm['concentration'].apply(pd.Series)
        del dfm['concentration']
        for col in concentration_columns:
            dfm[col] = [round(float(s), 2) for s in dfm[col].tolist()]
    else:
        dfm['concentration'] = [round(float(s), 2)
                                for s in dfm.concentration.tolist()]
    return dfm
