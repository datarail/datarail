import requests
import pandas as pd
import re


def get_ccle_mutations(cell_line='A549_LUNG'):
    """ Returns table of information on mutations in each cell line.
    Details include the Gene, chromsome, start and end poisition,
    and type of mutation

    Parameters
    ----------
    cell_line : str
      name of cell line and tissue (example: A549_LUNG)

    Returns
    -------
    df : pandas dataframe
       table of information on mutations found in input cell line
    """

    r = requests.get('https://portals.broadinstitute.org/ccle/'
                     'mutation_search/%s' % cell_line)
    df = pd.DataFrame(r.json()['data'])
    df.columns = ['Hugo_Symbol', 'Entrez_Gene_Id', 'NCBI_Build',
                  'Chromosome', 'Start_position', 'End_position',
                  'Variant_Classification', 'Varian_Type',
                  'Reference_Allele', 'Tumor_Seq_allele',
                  'dbSNP ID', 'Tumor Sample Barcode',
                  'Annotation_Transcript', 'Protein_Change',
                  'CCLE WES alt:ref', 'CCLE WGS alt:ref',
                  'CCLE RNAseq alt:ref', 'Sanger WES alt:ref ',
                  'CCLE HC1650 alt:ref', 'CCLE Raindance alt:ref',
                  'ExAC Allelic Fraction', 'TCGA Hotspot Count']
    df['Hugo_Symbol'] = [re.search(">(.*?)</a", hs).group(1)
                         for hs in df.Hugo_Symbol.tolist()]
    df['Entrez_Gene_Id'] = [re.search(">(.*?)</a", hs).group(1)
                            for hs in df.Entrez_Gene_Id.tolist()]
    return df
