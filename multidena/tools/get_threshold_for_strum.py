import pickle
import numpy
from strum import strum
from multidena.lib.common import read_seqs_with_complement, read_strum, get_threshold


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


def get_threshold_for_strum(fasta_path, strum_path, path_out):
    peaks = read_seqs_with_complement(fasta_path)
    strum_model = read_strum(strum_path)
    length = strum_model.k
    scores, number_of_sites = calculate_scores_strum_thresholds(peaks, strum_model, length)
    get_threshold(scores, number_of_sites, path_out)
    return(0)
