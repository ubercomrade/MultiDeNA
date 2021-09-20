import random
import shutil
import os
import subprocess
import bisect
import pickle
from shutil import copyfile
from operator import itemgetter
import numpy as np
from strum import strum
from lib.common import read_peaks, creat_background, \
write_fasta, complement, \
calculate_particial_auc, \
write_auc, calculate_merged_roc, \
calculate_short_roc, write_roc, \
calculate_fprs, \
make_pcm, make_pfm, write_meme, write_strum
from lib.speedup import creat_table_bootstrap


def strum_de_novo(fasta_path, model_length, cpu_count):
    fasta = open(fasta_path)
    strum_model = strum.StruM(mode='proteingroove', n_process=cpu_count)
    strum_model.train_EM(fasta, fasta=True, k=model_length,
                   lim=10**-5, n_init=5)
    fasta.close()
    return(strum_model)


def false_scores_strum(peaks, strum_model):
    false_scores = np.array([])
    peak = ''.join(peaks)
    complement_peak = complement(peak)
    false_scores = np.concatenate([false_scores,strum_model.score_seq(peak)])
    false_scores = np.concatenate([false_scores,strum_model.score_seq(complement_peak)])
    return(false_scores.tolist())


def true_scores_strum(peaks, strum_model, length):
    true_scores = []
    sites = []
    for peak in peaks:
        complement_peak = complement(peak)
        scores_1 = strum_model.score_seq(peak)
        scores_2 = strum_model.score_seq(complement_peak)
        index_1 = np.argmax(scores_1)
        index_2 = np.argmax(scores_2)
        if scores_1[index_1] > scores_2[index_2]:
            best = scores_1[index_1]
            site = peak[index_1:index_1+length]
        else:
            best = scores_2[index_2]
            site = complement_peak[index_2:index_2+length]
        true_scores.append(best)
        sites.append(site)
    return(true_scores, sites)


def fpr_at_tpr(true_scores, false_scores, tpr):
    true_scores.sort(reverse=True)
    false_scores.sort(reverse=True)
    false_length = len(false_scores)
    true_length = len(true_scores)
    score = true_scores[round(true_length * tpr) - 1]
    actual_tpr = sum([1 if true_score >= score else 0 for true_score in true_scores]) / true_length
    fpr = sum([1 if false_score >= score else 0 for false_score in false_scores]) / false_length
    return(fpr)


def write_sites(output, tag, sites):
    with open(output + '/' + tag + '.fasta', 'w') as file:
        for index, site in enumerate(sites):
            file.write(site + '\n')
    return(0)


def learn_optimized_strum(peaks_path, backgroud_path, counter, tmp_dir, output_auc, cpu_count, pfpr):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    if not os.path.isdir(output_auc):
        os.mkdir(output_auc)
    if os.path.exists(output_auc + '/auc.txt'):
        os.remove(output_auc + '/auc.txt')
    peaks = read_peaks(peaks_path)
    #for length in range(12, 41, 4):
    for length in range(10, 31, 2):
        true_scores = []
        false_scores = []
        sites = []
        for step in ['odd', 'even']:
            if step == 'odd':
                train_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 != 0]
                test_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 == 0]
            else:
                train_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 == 0]
                test_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 != 0]                
            train_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 != 0]
            test_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 == 0]
            write_fasta(train_peaks, tmp_dir + '/train.fasta')
            if os.path.isfile(backgroud_path):
                shuffled_peaks = read_peaks(backgroud_path)
            else:
                shuffled_peaks = creat_background(test_peaks, length, counter)
            strum_model = strum_de_novo(tmp_dir + '/train.fasta', length, cpu_count)
            for true_score, site in zip(*true_scores_strum(test_peaks, strum_model, length)):
                true_scores.append(true_score)
                sites.append(site)
            for false_score in false_scores_strum(shuffled_peaks, strum_model):
                false_scores.append(false_score)
        fprs = calculate_fprs(true_scores, false_scores)
        roc = calculate_short_roc(fprs, step=1)
        merged_roc = calculate_merged_roc(fprs)
        auc = calculate_particial_auc(merged_roc['TPR'], merged_roc['FPR'], pfpr)
        print("Length {};".format(length), "pAUC at {0} = {1};".format(pfpr, auc))
        write_auc(output_auc + '/auc.txt', auc, length)
        write_roc(output_auc + "/training_bootstrap_merged_{0}.txt".format(length), merged_roc)
        write_roc(output_auc + "/training_bootstrap_{0}.txt".format(length), roc)
        tag = 'strum_model_{0}'.format(length)
        sites = [i for i in sites if not 'N' in i]
        pcm = make_pcm(sites)
        pfm = make_pfm(pcm)
        nsites = len(sites)
        background = {'A': 0.25,
                     'C': 0.25,
                     'G': 0.25,
                     'T': 0.25}
        write_strum(strum_model, output_auc + '/{}.pickle'.format(tag))
        write_meme(output_auc, tag, pfm, background, nsites)
    shutil.rmtree(tmp_dir)
    return(0)


def choose_best_model(output_auc):
    auc = []
    with open(output_auc + '/auc.txt') as file:
        for line in file:
            auc.append(tuple(map(float, line.strip().split())))
        file.close()
    auc.sort(key=itemgetter(1))
    length = int(auc[-1][0])
    return(length)


def de_novo_with_oprimization_strum(peaks_path, backgroud_path, tmp_dir, output_dir, output_auc, cpu_count, pfpr):
    counter = 5000000
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    if not os.path.isdir(output_auc):
        os.mkdir(output_auc)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    learn_optimized_strum(peaks_path, backgroud_path, counter, tmp_dir, output_auc, cpu_count, pfpr)
    length = choose_best_model(output_auc)
    copyfile(output_auc + '/training_bootstrap_{}.txt'.format(length), 
             output_dir + '/bootstrap.txt')
    copyfile(output_auc + '/training_bootstrap_merged_{}.txt'.format(length), 
             output_dir + '/bootstrap_merged.txt')
    strum_model = strum_de_novo(peaks_path, length, cpu_count)
    peaks = read_peaks(peaks_path)
    true_scores, sites = true_scores_strum(peaks, strum_model, length)
    sites = [i for i in sites if not 'N' in i]
    pcm = make_pcm(sites)
    pfm = make_pfm(pcm)
    nsites = len(sites)
    background = {'A': 0.25,
                 'C': 0.25,
                 'G': 0.25,
                 'T': 0.25}
    tag = 'strum_model'
    write_meme(output_dir, tag, pfm, background, nsites)
    write_strum(strum_model, output_dir + '/{}.pickle'.format(tag))
    return(length)
