import numpy as np
import random


def creat_random_sample(test_sample, size_of_random_sample):
    container = []
    test_sample = [list(i) for i in test_sample]
    for index in range(len(test_sample) * size_of_random_sample):
        container.append(random.choice(test_sample))
        np.random.shuffle(container[-1])
    random_sites = []
    for site in container:
        random_sites.append(''.join(site))
    return(random_sites)


def creat_table_bootstrap(true_scores, false_scores):
    cdef list table = []
    true_scores.sort(reverse=True)
    false_scores.sort(reverse=True)
    false_length = len(false_scores)
    true_length = len(true_scores)
    for tpr in [round(i * 0.01, 2) for i in range(5,105, 5)]:
        score = true_scores[round(true_length * tpr) - 1]
        actual_tpr = sum([1 if true_score >= score else 0 for true_score in true_scores]) / true_length
        fpr = sum([1 if false_score >= score else 0 for false_score in false_scores]) / false_length
        table.append({'Scores': score, 'TPR': tpr, 'ACTUAL_TPR': actual_tpr, 'FPR': fpr})
    return(table)


def score_pwm(str seq, dict pwm):
    cdef float score = 0
    cdef int position = 0 
    cdef str letter
    for letter in seq:
        score += pwm[letter][position]
        position += 1
    return(score)


def calculate_scores_pwm_thresholds(list peaks, dict pwm, int length_of_site, float threshold):
    cdef str site
    cdef list scores = []
    cdef int i, N
    cdef int number_of_sites = 0
    append = scores.append
    for peak in peaks:
        N = len(peak) - length_of_site + 1
        for i in range(N):
            site = peak[i:length_of_site + i]
            if 'N' in site:
                continue
            number_of_sites += 1
            score = score_pwm(site, pwm)
            if score >= threshold:
                append(score)
    return(scores, number_of_sites)


def calculate_scores_pwm_bootstrap(list sites, dict pwm):
    cdef str site
    cdef list scores = []
    for site in sites:
        scores.append(score_pwm(site, pwm))
    return(scores)


def score_bamm(str site, dict bamm, int order, int length_of_site):
    cdef float score = 0
    cdef int index
    for index in range(order):
        score += bamm[site[0:index + 1]][index]
    for index in range(length_of_site - order):
        score += bamm[site[index:index+order + 1]][index + order]
    return(score)


def calculate_scores_bamm_bootstrap(list sites, dict bamm, int order, int length_of_site):
    cdef list scores = []
    append = scores.append
    for site in sites:
        append(score_bamm(site, bamm, order, length_of_site))
    return(scores)


def calculate_scores_bamm_thresholds(list peaks, dict bamm, int order, int length_of_site, float threshold):
    cdef list scores = []
    cdef int i, number_of_sites = 0
    cdef str peak
    append = scores.append
    for peak in peaks:
        for i in range(len(peak) - length_of_site + 1):
            site = peak[i:length_of_site + i]
            if 'N' in site:
                continue
            number_of_sites += 1
            score = score_bamm(site, bamm, order, length_of_site)
            if score >= threshold:
                append(score)
    return(scores, number_of_sites)