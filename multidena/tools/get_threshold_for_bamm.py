import itertools
from multidena.lib.common import read_seqs_with_complement, read_bamm, get_threshold
from multidena.lib.speedup import calculate_scores_bamm_thresholds


def min_score_bamm(bamm, order, length_of_site):
    scores = []
    min_score = 0
    for index in range(order):
        k_mers = itertools.product('ACGT', repeat=index + 1)
        scores.append([])
        for k in k_mers:
            k = ''.join(k)
            scores[-1].append(bamm[k][index])
    for index in range(length_of_site - order):
        k_mers = itertools.product('ACGT', repeat=order + 1)
        scores.append([])
        for k in k_mers:
            k = ''.join(k)
            scores[-1].append(bamm[k][index])
    for s in scores:
        min_score += min(s)
    return(min_score)


def max_score_bamm(bamm, order, length_of_site):
    scores = []
    max_score = 0
    for index in range(order):
        k_mers = itertools.product('ACGT', repeat=index + 1)
        scores.append([])
        for k in k_mers:
            k = ''.join(k)
            scores[-1].append(bamm[k][index])
    for index in range(length_of_site - order):
        k_mers = itertools.product('ACGT', repeat=order + 1)
        scores.append([])
        for k in k_mers:
            k = ''.join(k)
            scores[-1].append(bamm[k][index])
    for s in scores:
        max_score += max(s)
    return(max_score)


def to_score(norm_value, bamm, order, length_of_site):
    max_s = max_score_bamm(bamm, order, length_of_site)
    min_s = min_score_bamm(bamm, order, length_of_site)
    score = norm_value * (max_s - min_s) + min_s
    return(score)


def get_threshold_for_bamm(fasta_path, bamm_path, bg_path, path_out):
    peaks = read_seqs_with_complement(fasta_path)
    bamm, order = read_bamm(bamm_path, bg_path)
    length_of_site = len(bamm['A'])
    threshold = to_score(0.6, bamm, order, length_of_site)
    scores, number_of_sites = calculate_scores_bamm_thresholds(peaks, bamm, order, length_of_site, threshold)
    get_threshold(scores, number_of_sites, path_out)
    return(0)
