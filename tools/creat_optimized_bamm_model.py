import subprocess
import os
import shutil
import random
from lib.common import read_peaks, write_fasta, read_bamm, \
write_table_bootstrap, creat_background, \
score_bamm, complement, make_pcm, make_pfm, write_meme
from lib.speedup import creat_table_bootstrap


def create_bamm_model(peaks_path, directory, order, meme, extend, basename):
    args = ['BaMMmotif', directory, peaks_path, '--PWMFile', meme,
            '--CGS', '--order', str(order), '--Order', str(order),
           '--extend', str(extend), '--basename', str(basename)]
    r = subprocess.run(args, capture_output=True)
    bamm_path = directory + '/{}_motif_1.ihbcp'.format(basename)
    bg_path = directory + '/{}.hbcp'.format(basename)
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


def fpr_at_tpr(true_scores, false_scores, tpr):
    true_scores.sort(reverse=True)
    false_scores.sort(reverse=True)
    false_length = len(false_scores)
    true_length = len(true_scores)
    score = true_scores[round(true_length * tpr) - 1]
    actual_tpr = sum([1 if true_score >= score else 0 for true_score in true_scores]) / true_length
    fpr = sum([1 if false_score >= score else 0 for false_score in false_scores]) / false_length
    return(fpr)


def learn_optimized_bamm(peaks_path, counter, order, length, meme, tmp_dir):
    tpr = 0.3
    true_scores = []
    false_scores = []
    peaks = read_peaks(peaks_path)
    shuffled_peaks = creat_background(peaks, length, counter)
    bamm, order = create_bamm_model(peaks_path, tmp_dir, order, meme, 0, '0')
    for true_score in true_scores_bamm(peaks, bamm, order, length):
        true_scores.append(true_score)
    for false_score in false_scores_bamm(shuffled_peaks, bamm, order, length):
        false_scores.append(false_score)
    table = creat_table_bootstrap(true_scores, false_scores)
    fpr_current = fpr_at_tpr(true_scores, false_scores, tpr)
    index = 0
    print(length, fpr_current, order)
    for extend in range(1, 20):
        true_scores = []
        false_scores = []
        shuffled_peaks = creat_background(peaks, length + extend * 2, counter)
        bamm, order = create_bamm_model(peaks_path, tmp_dir, order, meme, extend, extend)
        for true_score in true_scores_bamm(peaks, bamm, order, length + extend):
            true_scores.append(true_score)
        for false_score in false_scores_bamm(shuffled_peaks, bamm, order, length + extend):
            false_scores.append(false_score)
        table = creat_table_bootstrap(true_scores, false_scores)
        fpr_new = fpr_at_tpr(true_scores, false_scores, tpr)
        if fpr_new < fpr_current: #and 1 - fpr_new/fprs_current < 0.5:
            fpr_current = fpr_new
            index += 1
            print(length + extend * 2, fpr_current, order)
        elif fpr_new > fpr_current and order <= 2:
            print(length + extend * 2, fpr_new, order)
            order += 1
        else:
            print(length + extend * 2, fpr_new, order)
            break
    return(order, index)


def de_novo_with_oprimization_bamm(peaks_path, length, meme, tmp_dir, output_dir):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    counter = 4000000
    order = 1
    bamm_order, model_index = learn_optimized_bamm(peaks_path, counter, order, length, meme, tmp_dir)
    shutil.copy(tmp_dir + '/{}_motif_1.ihbcp'.format(model_index),
           output_dir + '/bamm_model.ihbcp')
    shutil.copy(tmp_dir + '/{}.hbcp'.format(model_index),
           output_dir + '/bamm.hbcp')
    shutil.rmtree(tmp_dir)
    return(0)