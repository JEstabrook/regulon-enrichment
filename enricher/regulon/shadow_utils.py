from sklearn.feature_selection import f_regression, mutual_info_regression
from sklearn.mixture import BayesianGaussianMixture as GMM
from scipy.stats import spearmanr, pearsonr
from copy import copy
import warnings

warnings.simplefilter("ignore", UserWarning)
import pandas as pd
import dill as pickle
import functools
import os
import scipy.stats as st
import numpy as np
import random
import itertools

pd.options.mode.chained_assignment = None


def jaccard_sim(regulon, target1, target2):
    """ Calculate jaccard similarity coefficient on regulator pair space

    Args:
        regulon: obj subsetted regulon space
        target1: str First regulator in pair of regulators
        target2: str Second regulator in pair of regulators

    """
    a = set(regulon.loc[target1].Target)
    b = set(regulon.loc[target2].Target)
    u = a.union(b)
    c = a.intersection(b)
    # universe_overlap = float(len(c)) / (len(a) + len(b) - len(c))
    universe_overlap = float(len(c)) / (len(u))
    regul1_overlap = float(len(c)) / (len(a))
    regul2_overlap = float(len(c)) / (len(b))

    return universe_overlap, regul1_overlap, regul2_overlap, len(c), len(u), len(a), len(b)


def correlate_weights(regulon, target1, target2):
    """ Correlate regulon weights for pairs of regulators that overlap

    Args:
        regulon: obj subsetted regulon space
        target1: str First regulator in pair of regulators
        target2: str Second regulator in pair of regulators

    Returns:
        moa_corr: obj Dataframe containing correlation values of Mode of Activity
    """
    regul1 = regulon.loc[target1].set_index('Target')
    regul2 = regulon.loc[target2].set_index('Target')
    idx = regul1.index.intersection(regul2.index)
    sub_regul1 = regul1.reindex(idx)
    sub_regul2 = regul2.reindex(idx)
    moa_corr = sub_regul1.corrwith(sub_regul2)['MoA']
    return moa_corr


def calculate_overlap(enr, regulators = ['CEBPA', 'STAT5A', 'STAT5B', 'FOS', 'JUN', 'MYC', 'RUNX1'], n = 5):
    """ Calculates the percentage overalp of regualtor pairs

    Args:
        enr: obj enrichment object
        regulators (list): list regulators to perform shadow enrichment analysis on
        n: int number of shadow enrichment pairs to perform shadow analysis on

    Returns:

    """
    regulome = enr.regulon_weights.index.unique()
    permutations = list(itertools.product(regulators, regulome))
    res = ((*t, *jaccard_sim(enr.regulon_weights, *t), correlate_weights(enr.regulon_weights, *t)) for t in
           permutations)
    overlap = pd.DataFrame(res, columns = ['Regulator1', 'Regulator2', 'Universe', 'percent1', 'percent2', 'overlap',
                                           'edgecount', 'regulon_a_edgecount', 'regulon_b_edgecount',
                                           'correlation']).sort_values(
        ['Regulator1', 'correlation', 'Universe', 'percent1', 'percent2', 'overlap'],
        ascending = [False, False, False, False, False, False])
    shadow_regulators = overlap.groupby('Regulator1').head(n)
    selected_shadow_regulators = shadow_regulators[shadow_regulators.Universe != 1.0]

    return selected_shadow_regulators


def intersect_and_replace(enr, comb, seed):q
    """ Intersects regualtor pairs regulon space and replaces overlapping set with random edges

    Args:
        enr: obj enrichment object
        comb: tuple pair of regulators
        seed: int random seed

    Returns:
        boot_enrich: obj boostrap enrichment frame
    """
    random.seed(a = seed)
    subset_enr = copy(enr)
    # Determine features we can randomly sample
    universe = (set(enr.regulon.UpGene) | set(enr.regulon.DownGene)) & set(enr.expr.columns)
    # Subset the expression matrix
    subset_enr.expr = subset_enr.expr.reindex(universe, axis = 1)
    # Subset regulon by shadow enrichment pair
    regulon = subset_enr.regulon[subset_enr.regulon.UpGene.isin(comb)]
    # Adding subset below to see how this impacts enrichment
    regulon = regulon[regulon.UpGene.isin(enr.expr.columns) & regulon.DownGene.isin(enr.expr.columns)]
    # Identify features to be used in shadow enrichment for pair of regulators
    local_ = set(regulon.UpGene) | set(regulon.DownGene)
    # Remove the shared features from the universe feature set that will be randomly sampled
    diff = universe - local_
    regulon.reset_index(inplace = True, drop = True)
    # Identify the number of features that are represented in both regulons
    targs = (regulon.DownGene.value_counts() == 2)
    targs_list = targs[targs == True].index.tolist()
    try:
        # Try to randomly sample edges not represented in the two regulons
        random_targs = random.sample(diff, len(targs_list))
    except ValueError:
        # If the number of random edges to be sampled exceeds the total number of edges represented in [targs_list]
        # extend the number of edges to include other features
        sub_universe_ = (set(enr.regulon.UpGene) | set(enr.regulon.DownGene)) - local_ - diff
        targs_ = random.sample(sub_universe_, len(targs_list) - len(diff))
        print('Extending edge sample set by n={}'.format(len(targs_)))
        random_targs = list(targs_) + list(diff)

    map_dict = {k: v for k, v in zip(targs_list, random_targs)}
    regulon.DownGene.replace(map_dict, inplace = True)
    subset_enr.regulon = regulon
    subset_enr.assign_weights()
    quant_nes = regulon_enrichment.quantile_nes_score(subset_enr.regulon_weights, subset_enr.expr.T)
    nes_list = list(map(functools.partial(regulon_enrichment.score_enrichment, expr = subset_enr.expr,
                                          regulon = subset_enr.regulon_weights, quant_nes = quant_nes), comb))
    boot_enrich = pd.concat(nes_list, axis = 1)
    boot_enrich.columns = boot_enrich.columns + '_{}'.format(seed)

    return boot_enrich


def intersect_and_subset(enr, comb, seed):
    """ This function aims to perform an enrichment analysis for regulator TF1 and TF2 that co-regulate a number of
    genes. The enrichment score is calculated for TF1 and TF2 using the non-overlapping edges and then randomly
    resampling n regulated edges, where n is defined as the number of non-overlapping edges.

    Args:
        enr: enrichment obj
        comb: tuple pair of regulators
        seed: random seed int

    Returns:
        boot_enrich: enrichment obj
    """
    random.seed(a = seed)
    subset_enr = copy(enr)
    # Determine features we can randomly sample
    universe = (set(enr.regulon.UpGene) | set(enr.regulon.DownGene)) & set(enr.expr.columns)
    # Subset the expression matrix
    subset_enr.expr = subset_enr.expr.reindex(universe, axis = 1)
    # Subset regulon by shadow enrichment pair
    regulon = subset_enr.regulon[subset_enr.regulon.UpGene.isin(comb)]
    ### Adding subset below to see how this impacts enrichment
    regulon = regulon[regulon.UpGene.isin(enr.expr.columns) & regulon.DownGene.isin(enr.expr.columns)]
    # Identify features to be used in shadow enrichment for pair of regulators
    regulon.reset_index(inplace = True, drop = True)
    # Identify the number of features that are represented in both regulons
    targs = (regulon.DownGene.value_counts() == 2)
    targs_list = targs[targs == True].index.tolist()
    regul1 = regulon[regulon.UpGene.isin([comb[0]])]
    regul2 = regulon[regulon.UpGene.isin([comb[1]])]
    regul1_sub = regul1.sample(n = regul1.shape[0] - len(targs_list))
    regul2_sub = regul2.sample(n = regul2.shape[0] - len(targs_list))

    sub_regulon = pd.concat([regul1_sub, regul2_sub])

    subset_enr.regulon = sub_regulon
    subset_enr.assign_weights()
    quant_nes = regulon_enrichment.quantile_nes_score(subset_enr.regulon_weights, subset_enr.expr.T)
    nes_list = list(map(functools.partial(regulon_enrichment.score_enrichment, expr = subset_enr.expr,
                                          regulon = subset_enr.regulon_weights, quant_nes = quant_nes), comb))
    boot_enrich = pd.concat(nes_list, axis = 1)
    boot_enrich.columns = boot_enrich.columns + '_{}'.format(seed)

    return boot_enrich


def shadow_wrapper(combs, enr, nboot):
    """ Wraps intersect and replace function to generate dictionary that contains paired regulators and their
        recalculated enrichment signatures across n bootstraps
    Args:
        combs: list containing combination pairs of regulators
        enr: obj enrichment obj
        nboot: int number of boot straps

    Returns:
        stacked_dict: dict containing recalculated bootstrap enrichment scores for shadow pairs
    """
    shadow_dict = {comb: pd.concat([intersect_and_replace(enr, comb, n) for n in range(nboot)], axis = 1) for comb in
                   combs}
    reordered_dict = {
        comb: shadow_dict[comb].reindex(['{}_{}'.format(regulator, n) for regulator in comb for n in range(nboot)],
                                        axis = 1) for comb in combs}
    stacked_dict = {k: reordered_dict[k].stack().reset_index().rename(
        columns = {'level_0': 'Sample', 'level_1': 'boot', 0: 'Enrichment'}) for k in reordered_dict}
    return stacked_dict


def shadow_enrichment(enr, selected_shadow_regulators, nboot = 50):
    """ perform shadow enrichment on enrichment object using selected shadow regulators and n number of bootstraps

    Args:
        enr: obj enrichment object
        selected_shadow_regulators: obj dataframe outlining shadow pairs and regulon overlap
        nboot: int number of bootstraps to randomly resample regulon space

    Returns:
        shadow_dict: dict contains regulator pairs and their respective enrichment scores across bootstraps
    """
    if enr.total_enrichment is None:
        raise TypeError("`total_enrichment` must be assigned prior to shadow analysis")
    combs = [(x) for x, y in selected_shadow_regulators.groupby(['Regulator1', 'Regulator2'])]
    shadow_dict = shadow_wrapper(combs, enr, nboot)

    return shadow_dict
