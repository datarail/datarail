import pandas as pd


def combine_combos(df):
    dfc = df.groupby(['plate', 'well'], as_index=True).apply(
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
    dfi.columns = ['plate', 'well', 'agent', 'concentration']
    dfi['concentration'] = dfi['concentration'].fillna(0)
    dfm = combine_combos(dfi)
    del dfm['plate']
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
    return dfm
