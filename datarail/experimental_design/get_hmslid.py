import pandas as pd
# import re


def get_hmslid(drugs):
    ''' Given drug name, return hmslid'''

    url = 'http://lincs.hms.harvard.edu/db/sm/?search=&output_type=.csv'

    df = pd.read_csv(url)
    # alt_names = []
    # for i, names in enumerate(df['Alternative Names']):
    #     df['Alternative Names'][i] = re.split(';', names)
    #     alt_names += df['Alternative Names']

    hmslid = []

    for drug in drugs:
        if drug in df.Name.tolist():
            hmslid.append(df['HMS LINCS ID'].ix[df['Name'] == drug].values[0])
        elif drug in df['Alternative Names'].tolist():
            hmslid.append(df['HMS LINCS ID'].ix[
                df['Alternative Names'] == drug].values[0])
        else:
            hmslid.append('NA')

    drug_id = dict(zip(drugs, hmslid))

    return drug_id


