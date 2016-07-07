import os
import numpy as np
from pandas import read_csv, pivot_table, Series
from scipy.stats import hypergeom, norm
from statsmodels.sandbox.stats.multicomp import multipletests

__all__ = ['id_conversion', 'get_kinase_targets', 'get_example_data']


def id_conversion(lst, fr='uniprot', to='gene_name'):
    """
        Convert UniProt IDs to gene names (or the other way around),
            using a reduced version of the UniProt ID mapping table

        :param lst: List - List of identifiers, either UniProt IDs or Gene names or STRING-IDs
        :param fr: String - String indicating from which system the provided IDs are coming
        :param to: String - String indicating to which system should be matched;
                            For both parameters, combinations of 'uniprot', 'gene_name', or 'string_id' are allowed
        :return: List - List of converted identifiers
    """

    # Check for parameter consistency
    if fr == to:
        raise StandardError("Your function call doesn't make any sense: No need to match from %s to %s" % (fr, to))
    if len({fr, to}.intersection(['uniprot', 'gene_name', 'string_id'])) != 2:
        raise StandardError("Please provide only identifiers of supported systems!")

    # Load reduced id_mapping table from Uniprot
    id_mapping = read_csv(os.path.split(__file__)[0] + '/data/id_conversion.txt')

    # Map ids
    results = [id_mapping.ix[np.where(id_mapping[fr] == s)[0], to].values.tolist()[0]
               if len(id_mapping.ix[np.where(id_mapping[fr] == s)[0], to].values.tolist()) == 1 else np.nan
               for s in lst]
    return results


def get_kinase_targets(sources=None):
    """
        Retrieve kinase-substrate interactions from the OmniPath resource.
        :param sources: [String] - List of sources from which to import kinase-substrate interactions;
                                possible sources include:
                                ['ARN', 'CA1', 'dbPTM', 'DEPOD', 'HPRD', 'MIMP', 'Macrophage', 'NRF2ome',
                                'phosphoELM', 'PhosphoSite', 'SPIKE', 'SignaLink3', 'Signor', 'TRIP']

        :return: DataFrame - 1 for reported interaction (-1 for phosphatases), NaN for no interaction
                                    columns: kinase/phosphatase, rows: phospho-sites
    """

    # Check that only allowed sources are used
    if sources is None:
        sources = ['PhosphoSite']
    allowed_sources = ['ARN', 'CA1', 'dbPTM', 'DEPOD', 'HPRD', 'MIMP', 'Macrophage', 'NRF2ome',
                       'phosphoELM', 'PhosphoSite', 'SPIKE', 'SignaLink3', 'Signor', 'TRIP']
    if len(set(sources).difference(allowed_sources)) > 0 and sources != ['all']:
        raise StandardError('Please use only the sources integrated into pypath!')
    if sources == ['all']:
        sources = allowed_sources
    # Read in the data from OmniPath
    ptms_omnipath = read_csv(os.path.split(__file__)[0]+'/data/omnipath_ptms.txt', sep='\t')

    # Create unique identifier for each p-site
    ptms_omnipath['p_site'] = ptms_omnipath['UniProt_B'] + \
        '_' + \
        ptms_omnipath['Residue_letter'] + \
        ptms_omnipath['Residue_number'].astype(str)

    # Restrict the data on phosphorylation and dephosphorylation modification
    ptms_omnipath_phospho = ptms_omnipath.where(ptms_omnipath['PTM_type'].str.contains('phosphorylation')).dropna()

    # Add value column
    ptms_omnipath_phospho['value'] = 1

    if len(sources) != len(allowed_sources):
        # Include only interactions of specified sources
        ptms_omnipath_phospho = ptms_omnipath_phospho[[len(set(s.split(';')).intersection(sources)) >= 1
                                                       for s in ptms_omnipath_phospho['Databases']]]

    # Change value to negative value for dephosphorylation events
    ptms_omnipath_phospho.ix[ptms_omnipath_phospho['PTM_type'].str.startswith('de'), 'value'] = -1

    # Create pivot-table from pandas function
    adjacency_matrix = pivot_table(ptms_omnipath_phospho, values='value', index='p_site', columns='UniProt_A')

    # rename columns
    adjacency_matrix.columns = id_conversion(adjacency_matrix.columns.tolist())

    return adjacency_matrix


def __update_pypath_resource():
    """
        Update the file 'omnipath_ptms.txt',
        which is used to save the prior knowledge about kinase-substrate interactions from the Omnipath/pypath resource

        Use with care! Additionally, pypath must by installed on your machine.

        :return: None
    """
    # try importing pypath
    try:
        from pypath import main
    except:
        raise StandardError('The package pypath is not installed on your machine. '
                            'Please install pypath before you continue!')

    pa = main.PyPath()

    pa.init_network()

    pa.load_ptms()

    pa.export_ptms_tab(outfile=os.path.split(__file__)[0]+'/data/omnipath_ptms.txt')

    print 'pypath resource successfully updated!'


def get_example_data():
    """
        Read in the dataset from the publication from de Graaf et al. in MCP,
        and process it in order to use it for the other functions in kinact.

        :return: Tupel - (data_fc - fold change data organised as pandas DataFrame
                                    with different time points as columns and phosphosites as rows
                          data_p_value - p-values organised as data_fc, transformed as neg. logarithm with base 10
                          )
    """
    # Read data
    data_raw = read_csv(os.path.split(__file__)[0] + '/data/deGraaf_2014_jurkat.csv', sep=',', header=0)

    # Filter for those p-sites that were matched ambiguously
    data_reduced = data_raw[~data_raw['Proteins'].str.contains(';')]

    # Create identifier for each phosphorylation site, e.g. P06239_S59 for the Serine 59 in the protein Lck
    data_reduced.loc[:, 'ID'] = data_reduced['Proteins'] + '_' + data_reduced['Amino acid'] + \
        data_reduced['Positions within proteins']
    data_indexed = data_reduced.set_index('ID')

    # Extract only relevant columns
    data_relevant = data_indexed[[x for x in data_indexed if x.startswith('Average')]]

    # Rename columns
    data_relevant.columns = [x.split()[-1] for x in data_relevant]

    # Convert abundances into fold changes compared to control (0 minutes after stimulation)
    data_fc = data_relevant.sub(data_relevant['0min'], axis=0)
    data_fc.drop('0min', axis=1, inplace=True)

    # Also extract the p-values for the fold changes
    data_p_value = data_indexed[[x for x in data_indexed if x.startswith('p value') and x.endswith('vs0min')]]
    data_p_value.columns = [x.split('_')[-1].split('vs')[0] + 'min' for x in data_p_value]
    data_p_value = data_p_value.astype('float')  # Excel saved the p-values as strings, not as floating point numbers

    return data_fc, data_p_value
