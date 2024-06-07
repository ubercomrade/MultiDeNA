import csv
import itertools
import os
import sys
import math
import collections
from math import log
import re
import random
import bisect
import pickle
import operator
import pandas as pd
import numpy as np
import subprocess


# try:
#     from strum import strum
# except ModuleNotFoundError:
#     pass


def read_seqs(path):
    container = []
    letters = {'A', 'C', 'G', 'T'}
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                line = ''.join([l if l in letters else 'N' for l in line.strip().upper()])
                append(line.strip().upper())
    return(container)


def read_peaks(path):
    container = []
    letters = {'A', 'C', 'G', 'T'}
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                line = ''.join([l if l in letters else 'N' for l in line.strip().upper()])
                append(line.strip().upper())
    return(container)


def read_seqs_with_complement(path):
    container = []
    letters = {'A', 'C', 'G', 'T'}
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                line = ''.join([l if l in letters else 'N' for l in line.strip().upper()])
                append(line.strip().upper())
                append(complement(line.strip().upper()))
    return(container)


def complement(seq):
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


def read_fasta(path):
    fasta = list()
    letters = {'A', 'C', 'G', 'T'}
    with open(path, 'r') as file:
        for line in file:
            #print(line)
            if line.startswith('>'):
                line = line[1:].strip().split(':')
                record = dict()
                record['name'] = line[0]
                record['chromosome'] = line[2]
                coordinates_strand = line[3]

                start, end = re.findall(r'\d*-\d*', coordinates_strand)[0].split('-')
                record['start'] = start
                record['end'] = end

                strand = re.findall(r'\(.\)', coordinates_strand[:-3])
                if not strand == []:
                    record['strand'] = strand[0].strip('()')
                else:
                    record['strand'] = '+'
            else:
                line = ''.join([l if l in letters else 'N' for l in line.strip().upper()])
                record['seq'] = line.strip().upper()
                fasta.append(record)
    file.close()
    return(fasta)


def write_meme(output, tag, pfm, background, nsites):
    with open(output + '/' + tag + '.meme', 'w') as file:
        file.write('MEME version 4\n\nALPHABET= ACGT\n\nBackground letter frequencies\n')
        file.write('A {0} C {1} G {2} T {3}\n\n'.format(background['A'], background['C'],
                                                        background['G'], background['T']))
        file.write('MOTIF {0}\n'.format(tag))
        file.write(
            'letter-probability matrix: alength= 4 w= {0} nsites= {1}\n'.format(len(pfm['A']), nsites))
        for i in zip(pfm['A'], pfm['C'], pfm['G'], pfm['T']):
            file.write('{0:.8f}\t{1:.8f}\t{2:.8f}\t{3:.8f}\n'.format(i[0], i[1], i[2], i[3]))


def write_pwm(output, tag, pwm):
    with open(output + '/' + tag + '.pwm', 'w') as file:
        file.write('>{0}\n'.format(tag))
        for i in zip(pwm['A'], pwm['C'], pwm['G'], pwm['T']):
            file.write('{0}\t{1}\t{2}\t{3}\n'.format(i[0], i[1], i[2], i[3]))


def write_pfm(output, tag, pfm):
    with open(output + '/' + tag + '.pfm', 'w') as file:
        file.write('>{0}\n'.format(tag))
        for i in zip(pfm['A'], pfm['C'], pfm['G'], pfm['T']):
            file.write('{0:.9f}\t{1:.9f}\t{2:.9f}\t{3:.9f}\n'.format(i[0], i[1], i[2], i[3]))


def write_fasta(peaks, path):
    with open(path, 'w') as file:
        for index, p in enumerate(peaks):
            file.write('>{}\n'.format(index))
            file.write(p + '\n')
    return(0)


def creat_background(peaks, length_of_site, counter):
    shuffled_peaks = []
    number_of_sites = 0
    while counter > number_of_sites:
        peak = random.choice(peaks)
        shuffled_peak = ''.join(random.sample(peak, len(peak)))
        shuffled_peaks.append(shuffled_peak)
        number_of_sites += (len(''.join(shuffled_peak)) - length_of_site + 1) * 2
    return(shuffled_peaks)


# def write_table_bootstrap(path, data):
#     with open(path, 'w') as csvfile:
#         fieldnames = data[0].keys()
#         writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
#         writer.writeheader()
#         for line in data:
#             writer.writerow(line)
#     return(0)


def write_table(path, *args):
    with open(path, 'a') as file:
        file.write('\t'.join(list(map(str, args))[::-1])+'\n')
    return(0)


def write_roc(path, data):
    header = list(data.keys())
    with open(path, 'w') as file:
        file.write('\t'.join(header) + '\n')
        for i in zip(*data.values()):
            file.write('\t'.join((map(str,i))) + '\n')
    return(0)


def check_bootstrap(path):
    with open(path) as file:
        file.readline()
        for line in file:
            line = line.strip().split()
            tpr = float(line[1])
            fpr = float(line[3])
            if tpr == 0.5:
                break
    return(fpr)


def check_threshold_table(path):
    with open(path) as file:
        try:
            fpr = float(file.readline().strip().split()[1])
        except:
            fpr = -1
    return(fpr)


def write_scan(path, data):
    with open(path, 'w') as csvfile:
        fieldnames = ['chromosome', 'start', 'end', 'name', 'score', 'strand', 'site']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        #writer.writeheader()
        for line in data:
            writer.writerow(line)
    pass


# PWM MODEL
def make_pcm(motifs):
    matrix = {}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        matrix[i[0]] = []
    len_of_motif = len(motifs[0])
    for i in matrix.keys():
        matrix[i] = [0]*len_of_motif
    for i in range(len_of_motif):
        for l in motifs:
            matrix[l[i]][i] += 1
    return(matrix)


def make_pfm(pcm):
    number_of_sites = [0] * len(pcm['A'])
    for key in pcm.keys():
        for i in range(len(pcm[key])):
            number_of_sites[i] += pcm[key][i]
    pfm = dict()
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        pfm[i[0]] = []
    first_key = list(pcm.keys())[0]
    nuc_pseudo = 1/len(pcm.keys())
    for i in range(len(pcm[first_key])):
        for nuc in pcm.keys():
            pfm[nuc].append((pcm[nuc][i] + nuc_pseudo) / (number_of_sites[i] + 1))
    return(pfm)


def make_pwm(pfm):
    pwm = {}
    background = {'A': 0.25,
                 'C': 0.25,
                 'G': 0.25,
                 'T': 0.25}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        pwm[i[0]] = []
    first_key = list(pfm.keys())[0]
    for i in range(len(pfm[first_key])):
        for j in pfm.keys():
            pwm[j].append(math.log(pfm[j][i] / background[j]))
    return(pwm)


def sites_to_pwm(sites):
    pcm = make_pcm(sites)
    pfm = make_pfm(pcm)
    pwm = make_pwm(pfm)
    return(pwm)


def read_pwm(path):
    with open(path, 'r') as file:
        inf = file.readline()
        pwm = {'A': [], 'C': [], 'G': [], 'T': []}
        for line in file:
            line = line.strip().split('\t')
            for letter, value in zip(pwm.keys(), line):
                pwm[letter].append(float(value))
    file.close()
    return(pwm)


# diPWM MODEL

def make_dipcm(motifs):
    matrix = {}
    nucleotides = itertools.product('ACGT', repeat=2)
    for i in nucleotides:
        matrix["".join(i)] = []
    len_of_motif = len(motifs[0]) - 1
    for i in matrix.keys():
        matrix[i] = [0]*len_of_motif
    for i in range(len_of_motif):
        for l in motifs:
            matrix[l[i:i + 2]][i] += 1
    return(matrix)


def make_dipfm(dipcm):
    number_of_sites = [0] * len(dipcm['AA'])
    for key in dipcm.keys():
        for i in range(len(dipcm[key])):
            number_of_sites[i] += dipcm[key][i]
    dipfm = dict()
    nucleotides = itertools.product('ACGT', repeat=2)
    for i in nucleotides:
        dipfm["".join(i)] = []
    first_key = list(dipcm.keys())[0]
    nuc_pseudo = 1/len(dipcm.keys())
    for i in range(len(dipcm[first_key])):
        for dinuc in dipcm.keys():
            dipfm[dinuc].append((dipcm[dinuc][i] + nuc_pseudo) / (number_of_sites[i] + 1))
    return(dipfm)


def make_dipwm(dipfm):
    dipwm = {}
    background = 1 / 16
    nucleotides = itertools.product('ACGT', repeat=2)
    for i in nucleotides:
        dipwm["".join(i)] = []
    first_key = list(dipfm.keys())[0]
    for i in range(len(dipfm[first_key])):
        for j in dipfm.keys():
            dipwm[j].append(math.log(dipfm[j][i] / background))
    return(dipwm)


def sites_to_dipwm(sites):
    dipcm = make_dipcm(sites)
    dipfm = make_dipfm(dipcm)
    dipwm = make_dipwm(dipfm)
    return(dipwm)


# def read_dipwm(path):
#     with open(path, 'r') as file:
#         inf = file.readline()
#         pwm = {'A': [], 'C': [], 'G': [], 'T': []}
#         for line in file:
#             line = line.strip().split('\t')
#             for letter, value in zip(pwm.keys(), line):
#                 pwm[letter].append(float(value))
#     file.close()
#     return(pwm)


def read_dipwm(path):
    dipwm = {}
    nucleotides = itertools.product('ACGT', repeat=2)
    for i in nucleotides:
        dipwm["".join(i)] = []
    with open(path, 'r') as file:
        inf = file.readline()
        for line in file:
            line = line.strip().split('\t')
            for letter, value in zip(dipwm.keys(), line):
                dipwm[letter].append(float(value))
    file.close()
    return(dipwm)


def write_dipwm(output, tag, dipwm):
    container = []
    for i in dipwm.keys():
        container.append(dipwm[i])
    with open(output + '/' + tag + '.pwm', 'w') as file:
        file.write('>{0}\n'.format(tag))
        for i in zip(*container):
            file.write('\t'.join(map(str, i)) + '\n')


def write_dipfm(output, tag, dipfm):
    container = []
    for i in dipfm.keys():
        container.append(dipfm[i])
    with open(output + '/' + tag + '.pfm', 'w') as file:
        file.write('>{0}\n'.format(tag))
        for i in zip(*container):
            file.write('\t'.join(map(str, i)) + '\n')


def score_dipwm(seq, dipwm):
    score = 0
    position = 0
    length = len(seq) - 1
    for i in range(length):
        score += dipwm[seq[i:i+2]][i]
        position += 1
    return(score)


def calculate_scores_dipwm_thresholds(peaks, dipwm, length_of_site, threshold):
    scores = []
    number_of_sites = 0
    append = scores.append
    for peak in peaks:
        N = len(peak) - length_of_site + 1
        for i in range(N):
            site = peak[i:length_of_site + i]
            if 'N' in site:
                continue
            number_of_sites += 1
            score = score_dipwm(site, dipwm)
            if score >= threshold:
                append(score)
    return(scores, number_of_sites)


def calculate_scores_dipwm_bootstrap(sites, dipwm):
    scores = []
    for site in sites:
        scores.append(score_dipwm(site, dipwm))
    return(scores)


# BAMM MODEL
def parse_bamm_and_bg_from_file(bamm_file, bg_file):

    # Read BaMM file
    if os.path.isfile(bamm_file):
        motif_order = 0
        with open(bamm_file) as file:
            for line in file:
                if line[0] != '\n':
                    motif_order = motif_order + 1
                else:
                    break

        # count the motif length
        motif_length = int(sum(1 for line in open(bamm_file)) / (motif_order + 1))

        # read in bamm model
        model = {}
        for k in range(motif_order):
            model[k] = []

        with open(bamm_file) as file:
            for j in range(motif_length):
                for k in range(motif_order):
                    model[k].append([float(p) for p in file.readline().split()])
                file.readline()

    else:
        print('File {0} does not exist'.format(bg_file))
        sys.exit()

    # Read BG file
    bg = {}
    order = 0
    if os.path.isfile(bg_file):
        with open(bg_file) as bgmodel_file:
            line = bgmodel_file.readline()  # skip the first line for K
            line = bgmodel_file.readline()  # skip the second line for Alpha
            while order < motif_order:
                line = bgmodel_file.readline()
                bg_freq = [float(p) for p in line.split()]
                bg[order] = bg_freq
                order += 1
    else:
        print('File {0} does not exist'.format(bg_file))
        sys.exit()
    return(model, bg, order-1)


def make_k_mers(order):
    #  make list with possible k-mer based on bHMM model
    tmp = itertools.product('ACGT', repeat=order + 1)
    k_mer = []
    for i in tmp:
        k_mer.append(''.join(i))
    k_mer_dict = dict()
    index = 0
    for i in k_mer:
        k_mer_dict[i] = index
        index += 1
    return(k_mer_dict)


def make_log_odds_bamm(bamm, bg):
    log_odds_bamm = dict()
    for order in bamm.keys():
        log_odds_bamm[order] = [list(map(lambda x: log(x[0] / x[1], 2), zip(bamm_col, bg[order]))) for bamm_col in bamm[order]]
    return(log_odds_bamm)


def bamm_to_dict(log_odds_bamm, order, k_mers):
    bamm_dict = {}
    for k_mer in k_mers:
        bamm_dict[k_mer] = list()
    for i in range(len(log_odds_bamm[order])):
        for index, k_mer in enumerate(k_mers):
            bamm_dict[k_mer].append(log_odds_bamm[order][i][index])
    return(bamm_dict)


def read_bamm(bamm_path, bg_path):
    bamm, bg, order = parse_bamm_and_bg_from_file(bamm_path, bg_path)
    container = dict()
    log_odds_bamm = make_log_odds_bamm(bamm, bg)
    for i in range(order + 1):
        k_mers = make_k_mers(i)
        bamm_dict = bamm_to_dict(log_odds_bamm, i, k_mers)
        container.update(bamm_dict)
    return(container, order)


########################################################


def creat_random_sample(test_sample, size_of_random_sample):
    container = []
    test_sample = [list(i) for i in test_sample]
    for index in range(len(test_sample) * size_of_random_sample):
        container.append(random.choice(test_sample))
        random.shuffle(container[-1])
    random_sites = []
    for site in container:
        random_sites.append(''.join(site))
    return(random_sites)


def creat_table_bootstrap(true_scores, false_scores):
    table = []
    true_scores.sort(reverse=True)
    false_scores.sort(reverse=True)
    false_length = len(false_scores)
    true_length = len(true_scores)
    for tpr in [round(i * 0.01, 2) for i in range(1, 101, 1)]:
        score = true_scores[round(true_length * tpr) - 1]
        actual_tpr = sum([1 if true_score >= score else 0 for true_score in true_scores]) / true_length
        fpr = sum([1 if false_score >= score else 0 for false_score in false_scores]) / false_length
        table.append({'Scores': score, 'TPR': tpr, 'ACTUAL_TPR': actual_tpr, 'FPR': fpr})
    return(table)


def score_pwm(seq, pwm):
    score = 0
    position = 0
    for letter in seq:
        score += pwm[letter][position]
        position += 1
    return(score)


def calculate_scores_pwm_thresholds(peaks, pwm, length_of_site, threshold):
    scores = []
    number_of_sites = 0
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


def calculate_scores_pwm_bootstrap(sites, pwm):
    scores = []
    for site in sites:
        scores.append(score_pwm(site, pwm))
    return(scores)


def score_bamm(site, bamm, order, length_of_site):
    score = 0
    for index in range(order):
        score += bamm[site[0:index + 1]][index]
    for index in range(length_of_site - order):
        score += bamm[site[index:index+order + 1]][index + order]
    return(score)


def calculate_scores_bamm_bootstrap(sites, bamm, order, length_of_site):
    scores = []
    append = scores.append
    for site in sites:
        append(score_bamm(site, bamm, order, length_of_site))
    return(scores)


def calculate_scores_bamm_thresholds(peaks, bamm, order, length_of_site, threshold):
    scores = []
    number_of_sites = 0
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


# de-novo roc and auc


def calculate_short_roc(fprs, step=1):
    table = {'TPR': [0], 'FPR': [0]}#, 'SITES': [0]}
    current_number_of_sites = 0
    total_number_of_sites = len(fprs)
    for i in range(1, 100, step):
        position = round(total_number_of_sites * (i / 100))
        table['TPR'].append(i / 100)
        table['FPR'].append(fprs[position - 1])
        #table['SITES'].append(len(fprs[:position]))
    return(table)


def calculate_fool_roc(fprs):
    table = {'TPR': [0], 'FPR': [0]}#, 'SITES': [0]}
    current_number_of_sites = 1
    total_number_of_sites = len(fprs)
    for i in range(total_number_of_sites):
            table['TPR'].append(current_number_of_sites / total_number_of_sites)
            table['FPR'].append(fprs[i])
            #table['SITES'].append(current_number_of_sites)
            current_number_of_sites += 1
    return(table)


def calculate_merged_roc(fprs):
    table = {'TPR': [0], 'FPR': [0], 'SITES': [0]}
    current_number_of_sites = 0
    total_number_of_sites = len(fprs)
    for i in range(0, total_number_of_sites):
        if fprs[i] > fprs[i - 1]:
            table['TPR'].append(current_number_of_sites / total_number_of_sites)
            table['FPR'].append(fprs[i - 1])
            table['SITES'].append(current_number_of_sites)
        current_number_of_sites += 1
    table['TPR'].append(current_number_of_sites / total_number_of_sites)
    table['FPR'].append(fprs[i])
    table['SITES'].append(current_number_of_sites)
    return(table)


def calculate_roc_train(true_scores, false_scores):
    tprs = []
    fprs = []
    true_scores.sort()
    false_scores.sort()
    true_scores_uniq = list(set(true_scores))
    true_scores_uniq.sort(reverse=True)
    false_length = len(false_scores)
    true_length = len(true_scores)
    for score in true_scores_uniq:
        tpr = (true_length - bisect.bisect_right(true_scores, score)) / true_length
        fpr = (false_length - bisect.bisect_right(false_scores, score)) / false_length
        if fpr == 0:
            fpr = 0.5 / false_length
        tprs.append(tpr)
        fprs.append(fpr)
    return(tprs, fprs)


def calculate_fprs(true_scores, false_scores):
    fprs = []
    false_scores.sort()
    number_of_sites = len(false_scores)
    #true_scores_uniq = list(set(true_scores))
    #true_scores_uniq.sort(reverse=True)
    true_scores.sort(reverse=True)
    for score in true_scores:
        fpr = (number_of_sites - bisect.bisect_right(false_scores, score)) / number_of_sites
        if fpr == 0:
            fprs.append(0.5 / number_of_sites)
        else:
            fprs.append(fpr)
    return(fprs)


def calculate_prc(true_score, false_score):
    table = []
    prc = {'PRECISION': [0], 'RECALL': [0]}
    number_of_yes = len(true_score)
    yes = 0
    no = 0
    for t in true_score:
        table.append({'score': t, 'type': 1})
    for f in false_score:
        table.append({'score': f, 'type': 0})
    table.sort(key=operator.itemgetter('score'), reverse=True)
    score = table[0]['score']
    for line in table:
        if score != line['score']:
            score = line['score']
            precision = yes / (yes + no)
            recall = yes / number_of_yes
            prc['PRECISION'].append(precision)
            prc['RECALL'].append(recall)
        if line['type'] == 1:
            yes +=1
        else:
            no +=1
    precision = yes / (yes + no)
    recall = yes / number_of_yes
    prc['PRECISION'].append(precision)
    prc['RECALL'].append(recall)
    prc['PRECISION'][0] = prc['PRECISION'][1]
    return prc


def calculate_roc_bootstrap(fprs):
    table = [{'TPR': 0, 'FPR': 0, 'SITES': 0}]
    current_number_of_sites = 0
    total_number_of_sites = len(fprs)
    uniq_fprs = list(set(fprs))
    uniq_fprs.sort(reverse=True)
    counter = coollections.Counter(fprs)
    for val in enumerate(uniq_fprs):
        line = {}
        current_number_of_sites += counter[val]
        line['TPR'] = (current_number_of_sites / total_number_of_sites)
        line['FPR'] = (val)
        line['SITES'] = current_number_of_sites
        table.append(line)
    return(table)


def shorting_roc(roc):
    table = []
    for tpr in [round(i * 0.01, 2) for i in range(1, 100, 1)]:
        index = bisect.bisect_left(roc[0], tpr)
        actual_tpr = roc[0][index]
        fpr = roc[1][index]
        table.append({'TPR': tpr, 'ACTUAL_TPR': actual_tpr, 'FPR': fpr})
    return(table)


def write_auc(path, auc, length):
    with open(path, 'a') as file:
        file.write('{0}\t{1}\n'.format(length, auc))
    pass


def write_auc_with_order(path, auc, length, order):
    with open(path, 'a') as file:
        file.write('{0}\t{1}\t{2}\n'.format(order, length, auc))
    pass


def calculate_particial_auc(tprs, fprs, pfpr):
    auc = 0
    tpr_old = tprs[0]
    fpr_old = fprs[0]
    for tpr_new, fpr_new in zip(tprs[1:], fprs[1:]):
        if fpr_new >= pfpr:
            auc += (tpr_new + tpr_old) * ((pfpr - fpr_old) / 2)
            break
        auc += (tpr_new + tpr_old) * ((fpr_new - fpr_old) / 2)
        fpr_old = fpr_new
        tpr_old = tpr_new
    if fpr_new < pfpr:
        auc += (2* tpr_new) * ((pfpr - fpr_new) / 2)
    return(auc)


# Strum

def write_strum(strum_model, write_path):
    with open(write_path, 'wb') as f:
        pickle.dump(strum_model, f)
    pass


def read_strum(strum_path):
    with open(strum_path, 'rb') as f:
        strum_model = pickle.load(f)
    return(strum_model)


# Threshold table

def get_threshold(scores, number_of_sites, path_out):
    scores.sort(reverse=True) # big -> small
    with open(path_out, "w") as file:
        last_score = scores[0]
        last_fpr = 0
        for count, score in enumerate(scores[1:], 1):
            fpr = count/number_of_sites
            if score == last_score:
                continue
            elif count/number_of_sites > 0.001:
                file.write("{0}\t{1}\n".format(last_score, -np.log10(count/number_of_sites)))
                break
            elif score != last_score and fpr - last_fpr > 0.0000005:
                file.write("{0}\t{1}\n".format(last_score, -np.log10(count/number_of_sites))
                last_score = score
                last_fpr = fpr
    file.close()
    return(0)

# Annotation

#VG
# def split_by_chr(path, write_dir):
#     df = pd.read_csv(path, sep='\t', header=None)
#     df = df.sort_values(by=[0, 1])
#     chromosomes = set(df[0])
#     basename = os.path.basename(path).rsplit('.', 1)[0]
#     for chromosome in chromosomes:
#         chr_df = df[df[0] == chromosome]
#         chr_df.to_csv(f'{write_dir}/{basename}_{chromosome}.bed', sep='\t', header=None, index=False)
#     return chromosomes

#VG
# def gene_associated_with_motifs(path_scan, path_ann, tmp_dir, write_path):
#     basename_scan = os.path.basename(path_scan).rsplit('.', 1)[0]
#     basename_ann = os.path.basename(path_ann).rsplit('.', 1)[0]
#
#     chromosomes_scan = split_by_chr(path_scan, tmp_dir)
#     chromosomes_ann = split_by_chr(path_ann, tmp_dir)
#
#
#     for chromosome in chromosomes_ann:
#         if chromosome in chromosomes_scan:
#             args = ['/home/anton/Downloads/areajoin_v23/area_join_simple_one',
#                     f'{tmp_dir}/{basename_ann}_{chromosome}.bed',
#                     f'{tmp_dir}/{basename_scan}_{chromosome}.bed',
#                     f'{tmp_dir}/intersection_{chromosome}.bed',
#                     f'{tmp_dir}/aj.sta',
#                     '0']
#             r = subprocess.run(args, capture_output=True)
#
#
#     res = pd.DataFrame()
#     for chromosome in chromosomes_ann:
#         if chromosome in chromosomes_scan:
#             res = pd.concat([pd.read_csv(f'{tmp_dir}/intersection_{chromosome}.bed', sep='\t', header=None), res])
#     res = res.sort_values(by=[0, 1])
#     res = res[res[4] > 0]
#     res.to_csv(write_path, sep='\t', index=False, header=None)
#     return 0


# def gene_associated_with_motifs(path_scan, path_ann, tmp_dir, write_path):
#     args = ['bedtools', 'intersect',
#             '-a', path_ann,
#             '-b', path_scan,
#             '-wa', '-u',]
#     r = subprocess.run(args, capture_output=True)
#     res = pd.DataFrame([i.split() for i in r.stdout.decode().split('\n')])
#     res.to_csv(write_path, sep='\t', index=False, header=None)
#     return 0


# def gene_associated_with_motifs(path_scan, path_ann, write_path):
#     args = ['bedtools', 'intersect',
#             '-a', path_ann,
#             '-b', path_scan,
#             '-wa', '-u',]
#     r = subprocess.run(args, capture_output=True)
#     res = pd.DataFrame([i.split() for i in r.stdout.decode().strip().split('\n')])
#     if not res.empty:
#         genes = set(res[3])
#         with open(write_path, 'w') as file:
#             for gene in genes:
#                 file.write(f'{gene}\n')
#     else:
#         open(write_path, 'w').close()
#     return 0

def gene_associated_with_motifs(path_scan, path_ann, write_path):
    args = ['bedtools', 'intersect',
            '-a', path_ann,
            '-b', path_scan,
            '-wb',]
    r = subprocess.run(args, capture_output=True)
    res = pd.DataFrame([i.split() for i in r.stdout.decode().strip().split('\n')])
    if not res.empty:
        res.to_csv(write_path, sep='\t', index=False, header=None)
    else:
        open(write_path, 'w').close()
    #     genes = set(res[3])
    #     with open(write_path, 'w') as file:
    #         for gene in genes:
    #             file.write(f'{gene}\n')
    # else:
    #     open(write_path, 'w').close()
    return 0
