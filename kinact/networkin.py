from utils import get_kinase_targets, pivot_table, read_csv, os, np, Series, norm, multipletests


def prepare_networkin_files(phospho_sites, output_dir=os.getcwd() + '/networkin_files/'):
    """
        Prepare networkin prediction of up-stream kinases for a list of phospho-sites

        :param phospho_sites: List - List of phospho-site identifiers
        :param output_dir: Path - Directory in which to save the files that networkin needs

        :return: None - does not return anything, but saves the relevant files for the networkin analysis
    """
    # create output directory if it does not exist already
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # load UniProt/SwissProt sequences from the supplied data
    sequences = read_csv(os.path.split(__file__)[0] + '/data/sequences.tab', sep='\t', usecols=[0, 3], index_col=0)

    # open output files
    site_file = open(output_dir + 'site_file.txt', 'w+')
    fasta_file = open(output_dir + 'fasta_file.txt', 'w+')

    # keep track of which sequences have already been written in the fasta file
    added_sequences = list()

    for psite in phospho_sites:
        # enter each phospho-site in the sites-file in the required format for networkin
        site_file.write(psite.split('_')[0] + '\t' + psite.split('_')[1][1:] + '\t' + psite.split('_')[1][0] + '\n')

        # save the sequence of the complete protein into the fasta file, so that networkin can map the phospho-sites
        if not psite.split('_')[0] in added_sequences:
            if psite.split('_')[0] in sequences.index:
                fasta_file.write('>' + psite.split('_')[0] + '\n')
                fasta_file.write(sequences.ix[psite.split('_')[0]].values[0] + '\n')
                added_sequences.append(psite.split('_')[0])
            else:
                added_sequences.append(psite.split('_')[0])
    # Close files and print success message
    site_file.close()
    fasta_file.close()
    print 'Files for NetworKIN analysis successfully saved in %s' % output_dir


def get_kinase_targets_from_networkin(file_path, add_omnipath=True, score_cut_off=1):
    """
        Convert output of networkin to adjacency matrix

        :param file_path: Path - Path to the networkin output/result file
        :param add_omnipath: Boolean - Indicates whether to add the curated information from omnipath
        :param score_cut_off: Float - cut off for the networkin score

        :return: DataFrame - interactions between kinases/phosphatases and individual p-sites are assigned
                                    with their respective score from networkin
    """
    # read results file generate by networkin (release 3.0)
    nwkin_results = read_csv(file_path, sep='\t')

    # re-create unique phospho-site identifiers
    nwkin_results['p_site'] = nwkin_results['#Name'] + '_' + nwkin_results['Position']

    # restrict output to kinases and phosphatases
    nwkin_results = nwkin_results.where((nwkin_results['Tree'] == 'KIN') |
                                        (nwkin_results['Tree'] == 'PTP')).dropna(how='all')

    # Save phosphatases
    ptp = np.unique(nwkin_results['Kinase/Phosphatase/Phospho-binding domain description'].where(
        nwkin_results['Tree'] == 'PTP').dropna())

    # Create adjacency_matrix
    adjacency_matrix = pivot_table(data=nwkin_results,
                                   values='NetworKIN score',
                                   index='p_site',
                                   columns='Kinase/Phosphatase/Phospho-binding domain description')

    # Set all scores below the cut-off to NA
    adjacency_matrix = adjacency_matrix.where(adjacency_matrix > score_cut_off, np.nan)

    # add omnipath resources
    if add_omnipath:

        adjacency_matrix[adjacency_matrix > score_cut_off] = 1
        adjacency_matrix[adjacency_matrix < score_cut_off] = np.nan

        # Convert scores for phosphatases
        adjacency_matrix[ptp] = - adjacency_matrix[ptp]
        omnipath = get_kinase_targets(sources=['all'])

        adjacency_matrix = adjacency_matrix.reindex(columns=list(set(omnipath.columns.tolist()).union(
            adjacency_matrix.columns.tolist())))

        # Add information from curated data
        for p_site in list(set(omnipath.index).intersection(adjacency_matrix.index)):
            for kin in omnipath.loc[p_site, :].replace(0, np.nan).dropna().index:
                if omnipath.ix[p_site, kin] == -1:
                    adjacency_matrix.ix[p_site, kin] = -1
                else:
                    adjacency_matrix.ix[p_site, kin] = 1

        return adjacency_matrix
    else:

        # Convert scores for phosphatases
        adjacency_matrix[ptp] = - adjacency_matrix[ptp]

        return adjacency_matrix


def weighted_mean(data_fc, interactions, mP, delta, minimum_set_size=5):
    """
    Computes the kinase activity scores as mean/median of the fold changes in the substrate set.
    Statistical evaluation of the activity score is performed via a z-score.

    :param data_fc: Series - fold change data of a single condition as Pandas Series with phosphosites as Index
    :param interactions: DataFrame - Adjacency Matrix, for example created from kinact.get_kinase_targets()
    :param mP: Float - mean of the fold changes of the complete data set
    :param delta: Float - standard deviation of the complete data set
    :param minimum_set_size: Integer - minimum number of phosphosites present in substrate set of a kinase

    :return: Tupel - (scores - KSEA activity scores as Pandas Series,
                      p_value_adj - p-values of KSEA sccores, adjusted with Benjamini/Hochberg procedure
                      )
    """

    # Reduce interaction matrix to contain only phosphosites that are also detected in the data
    interactions_red = interactions.ix[list(set(data_fc.index.tolist()).intersection(interactions.index.tolist())), :]

    # Delete columns with less than minimum_set_size entries
    interactions_red = interactions_red.ix[:, interactions_red.replace(0, np.nan).notnull().sum() >= minimum_set_size]

    # Calculate mean of fold changes in susbtrate set, weighted with the networkin scores
    scores = Series({kinase: (data_fc*interactions_red[kinase]).sum()/float(interactions_red[kinase].replace(0, np.nan).dropna().sum())
                     for kinase in interactions_red})

    z_scores = Series({kinase: (scores[kinase]-mP) * np.sqrt(abs(interactions_red[kinase].replace(
        0, np.nan).sum())) * 1/delta for kinase in interactions_red})

    # Convert z-scores into p-values and adjust for multiple testing
    p_value = Series(norm.sf(z_scores), index=z_scores.index)
    # Adjust p-values for multiple testing using the Benjamini/Hochberg correction
    p_value_adj = Series(multipletests(p_value, alpha=0.05, method='fdr_bh')[1], index=p_value.index)

    return scores, p_value_adj