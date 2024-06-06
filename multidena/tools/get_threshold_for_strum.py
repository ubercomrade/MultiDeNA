import pickle
import numpy
from strum import strum
from multidena.lib.common import read_seqs_with_complement, read_strum


def calculate_scores_strum_thresholds(peaks, strum_model, length):
    merged_peaks = ''.join(peaks)
    real_indexes = []
    cumulative_length = 0
    for p in peaks:
        for i in range(cumulative_length, cumulative_length + len(p) - (length) + 1):
            real_indexes.append(i)
        cumulative_length += len(p)
    scores = strum_model.score_seq(merged_peaks)
    scores = scores[real_indexes]
    number_of_sites = len(scores)
    return(scores, number_of_sites)


def get_threshold(scores, number_of_sites, path_out):
    scores.sort(reverse=True) # big -> small
    with open(path_out, "w") as file:
        last_score = scores[0]
        for count, score in enumerate(scores[1:], 1):
            fpr = count/number_of_sites
            if score == last_score:
                continue
            elif count/number_of_sites > 0.001:
                file.write("{0}\t{1}\n".format(last_score, count/number_of_sites))
                break
            elif score != last_score and fpr - last_fpr > 0.0000005:
                file.write("{0}\t{1}\n".format(last_score, count/number_of_sites))
                last_score = score
                last_fpr = fpr
    file.close()
    return(0)



def get_threshold_for_strum(fasta_path, strum_path, path_out):
    peaks = read_seqs_with_complement(fasta_path)
    strum_model = read_strum(strum_path)
    length = strum_model.k
    scores, number_of_sites = calculate_scores_strum_thresholds(peaks, strum_model, length)
    get_threshold(scores, number_of_sites, path_out)
    return(0)
