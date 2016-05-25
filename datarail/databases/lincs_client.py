import pandas as pd


def get_drug_id(drugs):
    ''' Given drug name, return hmslid'''

    url = 'http://lincs.hms.harvard.edu/db/sm/?search=&output_type=.csv'

    df = pd.read_csv(url)

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


def get_cell_line_id(cell_lines):

    url = 'http://lincs.hms.harvard.edu/db/cells/?search=&output_type=.csv'

    df = pd.read_csv(url)

    hmslid = []

    for line in cell_lines:
        if line in df.Name.tolist():
            hmslid.append(df['HMS LINCS ID'].ix[df['Name'] == line].values[0])

        elif line in df['Alternative Names'].tolist():
            hmslid.append(df['HMS LINCS ID'].ix[
                df['Alternative Names'] == line].values[0])
        else:
            hmslid.append('NA')

    line_id = dict(zip(cell_lines, hmslid))

    return line_id
        

