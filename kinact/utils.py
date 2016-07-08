import os
import numpy as np
import warnings
from pandas import read_csv, pivot_table, Series, DataFrame
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


def get_kinase_targets(sources=None, organism='human'):
    """
        Retrieve kinase-substrate interactions from the OmniPath resource.
        :param sources: [String] - List of sources from which to import kinase-substrate interactions;
                                possible sources include:
                                ['ARN', 'CA1', 'dbPTM', 'DEPOD', 'HPRD', 'MIMP', 'Macrophage', 'NRF2ome',
                                'phosphoELM', 'PhosphoSite', 'SPIKE', 'SignaLink3', 'Signor', 'TRIP']
        :param organism: String - Organism of interest (currently: human, mouse, yeast)

        :return: DataFrame - 1 for reported interaction (-1 for phosphatases), NaN for no interaction
                                    columns: kinase/phosphatase, rows: phospho-sites
    """

    # Check if a supported organisms is supplied
    if organism != 'human' and organism != 'mouse' and organism != 'yeast':
        raise StandardError('No sequences available for the requested organism: %s' % organism)

    if organism == 'human':
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

    elif organism == 'yeast':
        # Check if sources had been supplied and raise warning if yes
        if sources is not None:
            warnings.warn('Sources parameter is only supported for human: value will be ignored!')
        # Read in data from phosphoGRID
        phospho_grid = read_csv(os.path.split(__file__)[0] + '/data/phospho_grid.txt', sep='\t', skiprows=32)

        # Create unique identifier
        phospho_grid['p_site'] = phospho_grid['A'] + '_' + phospho_grid['C']

        # Replace empty cells with np.nan
        phospho_grid.replace('-', np.nan, inplace=True)

        # Restrict to columns with known kinases (J) or phosphatases (N)
        kinases = phospho_grid[phospho_grid['K'].notnull()]
        phosphatases = phospho_grid[phospho_grid['O'].notnull()]

        # Extract interactions between all kinases and phosphatases
        kinases = [(site, kin, 1) for site, sources in kinases[['p_site', 'K']].values
                   for kin in sources.split('|') if kin != '']
        phosphatases = [(site, phos, -1) for site, sources in phosphatases[['p_site', 'O']].values
                        for phos in sources.split('|') if phos != '']
        interactions = kinases + phosphatases
        interactions = DataFrame(interactions, columns=['p_site', 'enzyme', 'value'])

        # Extract adjacency matrix
        adjacency_matrix = pivot_table(data=interactions, values='value',
                                       index='p_site', columns='enzyme', fill_value=np.nan)
        return adjacency_matrix

    elif organism == 'mouse':
        # Check if sources had been supplied and raise warning if yes
        if sources is not None:
            warnings.warn('Sources parameter is only supported for human: value will be ignored!')

        # Read in data from PhosphoSitePlus
        phosphosite = read_csv(os.path.split(__file__)[0] + '/data/PhosphoSitePlus.txt', sep='\t')
        phosphosite_mouse = phosphosite.loc[(phosphosite['KIN_ORGANISM'] == 'mouse') &
                                            (phosphosite['SUB_ORGANISM'] == 'mouse')]
        phosphosite_ksi = [(acc + '_' + rsd, kin)
                       for acc, rsd, kin in phosphosite_mouse[['SUB_ACC_ID', 'SUB_MOD_RSD', 'KINASE']].values]
        # and RegPhos
        reg_phos = read_csv(os.path.split(__file__)[0] + '/data/reg_phos_mouse.txt', sep='\t')
        reg_phos_red = reg_phos[reg_phos['catalytic kinase'].notnull()]
        reg_phos_ksi = [(acc + '_' + rsd + str(pos), kin.split('_group')[0])
                        if kin != 'autocatalysis'
                        else (acc + '_' + rsd + str(pos), sub_id.split('_MOUSE')[0])
                        for sub_id, acc, rsd, pos, kin
                        in reg_phos_red[['ID', 'AC', 'code', 'position', 'catalytic kinase']].values]
        # Combine interactions and construct adjacency matrix
        interactions = phosphosite_ksi + reg_phos_ksi
        interactions = DataFrame(interactions, columns=['p_site', 'kinase'])
        interactions['value'] = 1

        adjacency_matrix = pivot_table(data=interactions, index='p_site',
                                       columns='kinase', values='value', fill_value=np.nan)
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

    print 'Using pypath to update the prior-knowledge information about kinase-substrate interactions.'

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
