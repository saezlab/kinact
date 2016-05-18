#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# ######################################################################################################################
#
#  This file accompanies the chapter
#
#       'Phosphoproteomics-based profiling of kinase activities in cancer cells'
#
#  in the book
#
#       Methods of Molecular Biology: Cancer Systems Biology, Springer, 2016
#
#  in order to demonstrate the methodology and enable reproduction of the
#  kinase set enrichment analysis (KSEA) method, which was presented in
#
#       Casado et al.: Kinase-substrate enrichment analysis provides insights into the
#                      heterogeneity of signaling pathway activation in leukemia cells
#       Science signaling 268(6), 2013, doi: 10.1126/scisignal.2003573
#
#
#  Copyright (c): Queen Mary University of London (QMUL),
#                 Rheinisch-Westfälische Technische Hochschule (RWTH) Aachen
#
#  File author(s): Jakob Wirbel (jakob.wirbel@gmail.com)
#                  Emanuel Goncalves (emanuelgoncalves@gmail.com)
#
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
# ######################################################################################################################

# ######################################################################################################################
# Import required libraries
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, ttest_rel, norm
from statsmodels.sandbox.stats.multicomp import multipletests
# ######################################################################################################################

# ######################################################################################################################
# Functions for KSEA
# ######################################################################################################################


# ######################################################################################################################
# MEAN/MEDIAN: Function for KSEA calculates as mean of the fold changes of the substrate set
#
def ksea_mean(data_fc,              # the fold change data of a single condition, organised as a Pandas Series,
                                    # with the name of the p-sites as indices
              interactions,         # kinase-substrate relationships as adjacency matrix
                                    # (1 for interaction, 0 or NA otherwise)
              mP,                   # mean of the complete data set
              delta,                # standard deviation of the complete data set
              minimum_set_size=5,   # minimum number of p-sites present in substrate set of kinase, default is 5
              median=False          # Boolean indicating if mean or median should be used
              ):

    # Find intersection between kinase substrate sets and detected p-sites
    intersect = {kinase: list(set(interactions[kinase].replace(0, np.nan).dropna().index).intersection(data_fc.index))
                 for kinase in interactions}

    # Filter for kinases that have at least five target p-sites in the data
    intersect = {kinase: intersect[kinase] for kinase in intersect if len(intersect[kinase]) > minimum_set_size}

    # Global variables for the z-statistic
    mean_all = mP
    sd_all = delta

    if median:
        # Calculate the mean of the substrate set relative to the mean of the whole data set
        scores = pd.Series({kinase: np.median([data_fc.ix[intersect[kinase]]])
                            for kinase in intersect}).dropna()

        z_scores = pd.Series({kinase: abs((np.median([data_fc.ix[intersect[kinase]]]) - mean_all) *
                                          np.sqrt(len(intersect[kinase])) * 1 / sd_all)
                              for kinase in intersect}).dropna()
    else:
        # Calculate the mean of the substrate set relative to the mean of the whole data set
        scores = pd.Series({kinase: np.mean([data_fc.ix[intersect[kinase]]])
                            for kinase in intersect}).dropna()

        z_scores = pd.Series({kinase: abs((np.mean([data_fc.ix[intersect[kinase]]]) - mean_all) *
                                          np.sqrt(len(intersect[kinase])) * 1 / sd_all)
                              for kinase in intersect}).dropna()

    # Convert z-scores into p-values and adjust for multiple testing
    p_value = pd.Series(norm.sf(z_scores), index=z_scores.index)
    # Adjust p-values for multiple testing using the Benjamini/Hochberg correction
    p_value_adj = pd.Series(multipletests(p_value, alpha=0.05, method='fdr_bh')[1], index=p_value.index)

    return scores, p_value_adj


# ######################################################################################################################
# DELTA: Function for KSEA with number of significantly up-regulated minus significantly down-regulated p-sites
#
def ksea_delta(data_fc,                  # the fold change data of a single condition, organised as a Pandas Series,
                                         # with the name of the p-sites as indices
               p_values,                 # the p-values of an statistical test if the fold change is significant
                                         # compared to the control; it is assumed that the p-values have been
                                         # transformed with -log10(p-value)
               interactions,             # kinase-substrate relationships as adjacency matrix
                                         # (1 for interaction, 0 or NA otherwise)
               cut_off=-np.log10(0.05),  # cut-off for the p-value to define significant regulation, default is 0.05
               minimum_set_size=5,       # minimum number of p-sites present in substrate set of kinase, default is 5
               ):

    # Global parameters for the hyper-geometric test
    n_total = len(data_fc.dropna())
    n_reg = len(np.where(p_values.ix[interactions.index.tolist()].dropna() > cut_off)[0])

    # Find intersection between kinase substrate sets and detected p-sites
    intersect = {kinase: list(set(interactions[kinase].replace(0, np.nan).dropna().index).intersection(data_fc.index))
                 for kinase in interactions}

    # Filter for kinases that have at least five target p-sites in the data
    intersect = {kinase: intersect[kinase] for kinase in intersect if len(intersect[kinase]) > minimum_set_size}

    # Initialise hypergeometric distribution for each kinase
    hyper_geom_dist = {kinase: hypergeom(n_total, len(intersect[kinase]), n_reg) for kinase in intersect}

    # Calculate the number of p-sites in the substrate group that are significantly increased
    #                                                       minus the ones that are decreased
    scores = pd.Series({kinase: len(data_fc.ix[intersect[kinase]].where(
        (data_fc.ix[intersect[kinase]] > 0) & (p_values.ix[intersect[kinase]] > cut_off)).dropna()) -
                                len(data_fc.ix[intersect[kinase]].where(
        (data_fc.ix[intersect[kinase]] < 0) & (p_values.ix[intersect[kinase]] > cut_off)).dropna())
                        for kinase in intersect})

    # Calculate p-value of the score using the hyper-geometric distribution
    p_value = pd.Series({kinase: hyper_geom_dist[kinase].pmf(len(np.where(p_values.ix[intersect[kinase]] > cut_off)[0]))
                        if len(np.where(p_values.ix[intersect[kinase]] > cut_off)[0]) > 0 else 1
                        for kinase in intersect})
    # Benjamini-Hochberg correction of the p-values
    p_value_adj = pd.Series(multipletests(p_value, alpha=0.05, method='fdr_bh')[1], index=p_value.index)

    return scores, p_value_adj


# ######################################################################################################################
# MEAN - Alternative: Function for KSEA with mean of substrate set, considering only significantly changing p-sites
#
def ksea_mean_alt(data_fc,                  # the fold change data of a single condition, organised as a Pandas Series,
                                            # with the name of the p-sites as indices
                  p_values,                 # the p-values of an statistical test if the fold change is significant
                                            # compared to the control; it is assumed that the p-values have been
                                            # transformed with -log10(p-value)
                  interactions,             # kinase-substrate relationships as adjacency matrix
                                            # (1 for interaction, 0 or NA otherwise)
                  mP,                       # mean of the complete data set
                  delta,                    # standard deviation of the complete data set
                  cut_off=-np.log10(0.05),  # cut-off for the p-value to define significant regulation, default is 0.05
                  minimum_set_size=5,       # minimum number of p-sites present in substrate set of kinase, default is 5
                  median=False              # Boolean indicating if mean or median should be used
                  ):

    # Find intersection between kinase substrate sets and detected p-sites
    intersect = {kinase: list(set(interactions[kinase].replace(0, np.nan).dropna().index).intersection(data_fc.index))
                 for kinase in interactions}

    # Filter for kinases that have at least five target p-sites in the data
    intersect = {kinase: intersect[kinase] for kinase in intersect if len(intersect[kinase]) > minimum_set_size}

    reduced_substrate_set = {kinase: data_fc.ix[intersect[kinase]].where(
        p_values.ix[intersect[kinase]] > cut_off).dropna()
                             for kinase in intersect}

    # Global variables for the z-statistic
    mean_all = mP
    sd_all = delta

    if median:
        # Calculate the median of the substrate set considering only significantly changed p-sites
        scores = pd.Series({kinase: np.median(reduced_substrate_set[kinase]) for kinase in intersect}).dropna()

        z_scores = pd.Series(
            {kinase: abs((np.median(reduced_substrate_set[kinase]) - mean_all) *
                     np.sqrt(len(reduced_substrate_set[kinase])) * 1 / sd_all)
             for kinase in intersect}).dropna()
    else:
        # Calculate the mean of the substrate set considering only significantly changed p-sites
        scores = pd.Series({kinase: np.mean(reduced_substrate_set[kinase]) for kinase in intersect}).dropna()

        z_scores = pd.Series(
            {kinase: abs((np.mean(reduced_substrate_set[kinase]) - mean_all) *
                         np.sqrt(len(reduced_substrate_set[kinase])) * 1 / sd_all)
             for kinase in intersect}).dropna()

    # Convert z-scores into p-values and adjust for multiple testing
    p_value = pd.Series(norm.sf(z_scores), index=z_scores.index)
    # Adjust p-values for multiple testing using the Benjamini/Hochberg correction
    p_value_adj = pd.Series(multipletests(p_value, alpha=0.05, method='fdr_bh')[1], index=p_value.index)

    return scores, p_value_adj
