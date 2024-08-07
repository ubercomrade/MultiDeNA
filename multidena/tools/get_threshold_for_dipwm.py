from multidena.lib.common import read_seqs_with_complement, read_dipwm, get_threshold
from multidena.lib.speedup import calculate_scores_dipwm_thresholds


def to_score(norm_value, dipwm):
    min_s = min_score(dipwm)
    max_s = max_score(dipwm)
    score = norm_value * (max_s - min_s) + min_s
    return(score)


def to_norm(score, dipwm):
    min_s = min_score(dipwm)
    max_s = max_score(dipwm)
    norm_value = (score - min_s) / (max_s - min_s)
    return(norm_value)


def min_score(dipwm):
    value = int()
    keys = list(dipwm.keys())
    length_dipwm = len(dipwm[keys[0]])
    for i in range(length_dipwm):
        tmp = []
        for j in keys:
            tmp.append(dipwm[j][i])
        value += min(tmp)
    return(value)


def max_score(dipwm):
    value = int()
    keys = list(dipwm.keys())
    length_dipwm = len(dipwm[keys[0]])
    for i in range(length_dipwm):
        tmp = []
        for j in keys:
            tmp.append(dipwm[j][i])
        value += max(tmp)
    return(value)


def get_threshold_for_dipwm(fasta_path, dipwm_path, path_out):
    peaks = read_seqs_with_complement(fasta_path)
    dipwm = read_dipwm(dipwm_path)
    length_of_site = len(dipwm['AA']) + 1
    threshold = to_score(0.6, dipwm)
    scores, number_of_sites = calculate_scores_dipwm_thresholds(peaks, dipwm, length_of_site, threshold)
    get_threshold(scores, number_of_sites, path_out)
    return(0)
