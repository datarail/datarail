import requests
import pandas as pd
import re
import os


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


def get_hmslincs_annotations(cell_lines):
    """ Returns metadata annotations for cell lines available
    on LINCS DB
    Parameters
    ----------
    cell_lines: list of str
       Name of cell lines
    Return
    ------
    dfs: pandas dataframe
      dataframe of cell line annotations from HMS LINCS DB
    for given list of cell lines
    """
    dfcell = pd.read_csv('http://lincs.hms.harvard.edu/db/'
                         'cells/?search=&output_type=.csv')
    dfcell.index = dfcell.Name.tolist()
    relcols = ['HMS LINCS ID',  'Name',  'Alternative Names',
               'LINCS ID', 'Organ', 'Cell Type', 'Known Mutations',
               'Disease', 'Details of Disease',
               'Growth Properties', 'Recommended Culture Conditions']
    dfs = dfcell.loc[cell_lines, relcols]
    return dfs


def get_ccle_cnv(cell_lines, genes):
    """Returns copy number variations from CCLE for given list
    of cell lines and genes.
    Parameters
    ----------
    cell_lines : list of str
         cell_lines of interest.
         Naming convention should be consistent with CCLE.
    genes : list of str
         genes of interest (gene symbols).
    Returns
    -------
    dfs : pandas datafame
         data table of cnv for given list of cell lines and genes.
    """
    cla_path = os.path.dirname(os.path.abspath(__file__))
    cnv_file_path = os.path.join(cla_path,
                                 'CCLE_copynumber_byGene_2013-12-03.txt')
    dfcnv = pd.read_table(cnv_file_path,
                          index_col=0)
    dfs = dfcnv.loc[genes, cell_lines].copy()
    return dfs
