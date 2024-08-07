from multidena.lib.common import read_seqs_with_complement, read_pwm, get_threshold
from multidena.lib.speedup import calculate_scores_pwm_thresholds


def to_score(norm_value, pwm):
    min_s = min_score(pwm)
    max_s = max_score(pwm)
    score = norm_value * (max_s - min_s) + min_s
    return(score)


def to_norm(score, pwm):
    min_s = min_score(pwm)
    max_s = max_score(pwm)
    norm_value = (score - min_s) / (max_s - min_s)
    return(norm_value)


def min_score(pwm):
    value = int()
    keys = list(pwm.keys())
    length_pwm = len(pwm[keys[0]])
    for i in range(length_pwm):
        tmp = []
        for j in keys:
            tmp.append(pwm[j][i])
        value += min(tmp)
    return(value)


def max_score(pwm):
    value = int()
    keys = list(pwm.keys())
    length_pwm = len(pwm[keys[0]])
    for i in range(length_pwm):
        tmp = []
        for j in keys:
            tmp.append(pwm[j][i])
        value += max(tmp)
    return(value)


def get_threshold_for_pwm(fasta_path, pwm_path, path_out):
    peaks = read_seqs_with_complement(fasta_path)
    pwm = read_pwm(pwm_path)
    length_of_site = len(pwm['A'])
    threshold = to_score(0.6, pwm)
    scores, number_of_sites = calculate_scores_pwm_thresholds(peaks, pwm, length_of_site, threshold)
    get_threshold(scores, number_of_sites, path_out)
    return(0)
