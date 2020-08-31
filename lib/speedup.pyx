import random


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
    cdef int length = len(seq)
    for index in range(length):
        score += pwm[seq[index]][position]
        position += 1
    return(score)


def calculate_scores_pwm_thresholds(list peaks, dict pwm, int length_of_site, float threshold):
    cdef str site
    cdef list scores = []
    cdef int i, N
    cdef int number_of_sites = 0
    cdef int number_of_peaks = len(peaks)
    append = scores.append
    for index in range(number_of_peaks):
        peak = peaks[index]
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


def score_bamm(str site, dict bamm, int order, int length_of_site):
    cdef float score = 0
    cdef int index
    for index in range(order):
        score += bamm[site[0:index + 1]][index]
    for index in range(length_of_site - order):
        score += bamm[site[index:index+order + 1]][index + order]
    return(score)


def calculate_scores_bamm_thresholds(list peaks, dict bamm, int order, int length_of_site, float threshold):
    cdef list scores = []
    cdef int i, number_of_sites = 0
    cdef str peak
    cdef int number_of_peaks = len(peaks)
    append = scores.append
    for index in range(number_of_peaks):
        peak = peaks[index]
        for i in range(len(peak) - length_of_site + 1):
            site = peak[i:length_of_site + i]
            if 'N' in site:
                continue
            number_of_sites += 1
            score = score_bamm(site, bamm, order, length_of_site)
            if score >= threshold:
                append(score)
    return(scores, number_of_sites)