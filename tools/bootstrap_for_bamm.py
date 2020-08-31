import subprocess
import os
import shutil
import random
from lib.common import read_peaks, write_fasta, read_bamm, \
write_table_bootstrap, creat_background, score_bamm, complement
from lib.speedup import creat_table_bootstrap


def create_bamm_model(directory, order, meme):
    fasta_path = directory + '/train.fasta'
    args = ['BaMMmotif', directory, fasta_path, '--PWMFile', meme, '--EM', '--order', str(order), '--Order', str(order)]
    r = subprocess.run(args, capture_output=True)
    bamm_path = directory + '/train_motif_1.ihbcp'
    bg_path = directory + '/train.hbcp'
    bamm, order = read_bamm(bamm_path, bg_path)
    return(bamm, order)


def false_scores_bamm(peaks, bamm, order, length_of_site):
    false_scores = []
    append = false_scores.append
    for peak in peaks:
        complement_peak = complement(peak)
        full_peak = peak + 'N' * length_of_site + complement_peak
        n = len(full_peak) - length_of_site + 1
        for i in range(n):
            site = full_peak[i:length_of_site + i]
            if 'N' in site:
                continue
            score = score_bamm(site, bamm, order, length_of_site)
            false_scores.append(score)
    return(false_scores)


def true_scores_bamm(peaks, bamm, order, length_of_site):
    true_scores = []
    for peak in peaks:
        complement_peak = complement(peak)
        best = -1000000
        full_peak = peak + 'N' * length_of_site + complement_peak
        n = len(full_peak) - length_of_site + 1
        for i in range(n):
            site = full_peak[i:length_of_site + i]
            if 'N' in site:
                continue
            score = score_bamm(site, bamm, order, length_of_site)
            if score >= best:
                best = score
        true_scores.append(best)
    return(true_scores)


def bootstrap_bamm(peaks, length_of_site, counter, order, meme, tmp_dir):
    true_scores = []
    false_scores = []
    number_of_peaks = len(peaks)
    for i in range(10):
        train_peaks = random.choices(peaks, k=round(0.9 * number_of_peaks))
        test_peaks = [peak for peak in peaks if not peak in train_peaks]
        shuffled_peaks = creat_background(test_peaks, length_of_site, counter / 10)
        write_fasta(train_peaks, tmp_dir + '/train.fasta')
        bamm, order = create_bamm_model(tmp_dir, order, meme)
        for true_score in true_scores_bamm(test_peaks, bamm, order, length_of_site):
            true_scores.append(true_score)
        for false_score in false_scores_bamm(shuffled_peaks, bamm, order, length_of_site):
            false_scores.append(false_score)
    table = creat_table_bootstrap(true_scores, false_scores)
    return(table)


def bootstrap_for_bamm(peaks_path, results_path, length_of_site, meme, tmp_dir, counter = 5000000, order=2):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    peaks = read_peaks(peaks_path)
    table = bootstrap_bamm(peaks, length_of_site, counter, order, meme, tmp_dir)
    write_table_bootstrap(results_path, table)
    shutil.rmtree(tmp_dir)
    return(0)


