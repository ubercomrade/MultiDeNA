import csv
import itertools
import os
import sys
import math
import re
import random
import subprocess
import shutil
import glob
from operator import itemgetter
import argparse
from lib.common import calculate_roc, calculate_particial_auc


def read_peaks(path):
    container = []
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                append(line.strip().upper())
    return(container)


def read_inmode_bed(path):
    table = []
    with open(path) as file:
        for line in file:
            line = line.strip().split('\t')
            line[0] = int(line[0])
            line[1] = int(line[1])
            line[2] = int(line[2])
            line[4] = float(line[4])
            table.append(line)
    return(table)


def write_fasta(sites, tmp_dir, tag):
    with open('{0}/{1}.fa'.format(tmp_dir, tag), 'w') as file:
        for index, site in enumerate(sites):
            file.write('>{0}\n{1}\n'.format(index, site))
    return(0)


def true_scores_inmode(path_to_inmode, path_to_java, motif_length, tmp_dir, fasta_tag, model_tag):
    scores = []
    args = [path_to_java, '-Xmx16G', '-Xms1G', 
            '-jar',
            path_to_inmode, 'scan',
            'i={0}/{1}_inmode_model.xml'.format(tmp_dir, model_tag),
            'id={0}/{1}.fa'.format(tmp_dir, fasta_tag), 'f=1.0', 'outdir={}'.format(tmp_dir), 'bs=false']
    r = subprocess.run(args, capture_output=True)
    scores = []
    table = read_inmode_bed('{0}/{1}'.format(tmp_dir, "/Motif_hits_from_SequenceScan(1.0).BED"))
    table.sort(key=itemgetter(0, 4))
    last_index = 0
    for line in table:
        index = line[0]
        score = line[4]
        if last_index != index:
            scores.append(last_score)
        last_score = score
        last_index = index
    scores.append(score)
    scores = [math.log(float(i), 10) for i in scores]
    os.remove(tmp_dir + '/Binding_sites_from_SequenceScan(1.0).txt')
    os.remove(tmp_dir + '/Motif_hits_from_SequenceScan(1.0).BED')
    os.remove(tmp_dir + '/protocol_scan.txt')
    #os.remove(tmp_dir + '/{}.fa'.format(tag))
    return(scores)


def false_scores_inmode(path_to_inmode, path_to_java, motif_length, tmp_dir, fasta_tag, model_tag):
    scores = []
    args = [path_to_java, '-Xmx16G', '-Xms1G', 
            '-jar',
            path_to_inmode, 'scan',
            'i={0}/{1}_inmode_model.xml'.format(tmp_dir, model_tag),
            'id={0}/{1}.fa'.format(tmp_dir, fasta_tag), 'f=1.0', 'outdir={}'.format(tmp_dir), 'bs=false']
    r = subprocess.run(args, capture_output=True)
    with open('{0}/{1}'.format(tmp_dir, "/Motif_hits_from_SequenceScan(1.0).BED")) as file:
        for line in file:
            scores.append(math.log(float(line.split()[4]), 10))
    os.remove(tmp_dir + '/Binding_sites_from_SequenceScan(1.0).txt')
    os.remove(tmp_dir + '/Motif_hits_from_SequenceScan(1.0).BED')
    os.remove(tmp_dir + '/protocol_scan.txt')
    #os.remove(tmp_dir + '/{}.fa'.format(fasta_tag))
    return(scores)


def make_inmode(peaks_path, path_to_inmode, path_to_java, motif_length, order, tmp_dir, tag):
    args = [path_to_java, '-Xmx16G', '-Xms1G', '-jar', path_to_inmode,
    'denovo', 'i={}'.format(peaks_path), 'm={}'.format(motif_length), 'outdir={}'.format(tmp_dir),
    'mo={}'.format(order)]
    r = subprocess.run(args, capture_output=True)
    shutil.copy(tmp_dir + '/Learned_DeNovo({0},{1},2)_motif/XML_of_DeNovo({0},{1},2)_motif.xml'.format(motif_length,order),
                tmp_dir + '/{}_inmode_model.xml'.format(tag))
    os.remove(tmp_dir + '/protocol_denovo.txt')
    os.remove(tmp_dir + '/Logfile_of_DeNovo({0},{1},2)_stochastic_search.txt'.format(motif_length,order))
    os.remove(tmp_dir + '/Latent_variables_of_DeNovo({0},{1},2).txt'.format(motif_length,order))
    shutil.rmtree(tmp_dir + '/Learned_DeNovo({0},{1},2)_motif'.format(motif_length,order))
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


def fpr_at_tpr(true_scores, false_scores, tpr):
    true_scores.sort(reverse=True)
    false_scores.sort(reverse=True)
    false_length = len(false_scores)
    true_length = len(true_scores)
    score = true_scores[round(true_length * tpr) - 1]
    actual_tpr = sum([1 if true_score >= score else 0 for true_score in true_scores]) / true_length
    fpr = sum([1 if false_score >= score else 0 for false_score in false_scores]) / false_length
    return(fpr)


def complement(seq):
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


def learn_optimized_inmode(peaks_path, counter, order, length, 
    path_to_inmode, path_to_java, tmp_dir, tpr, pfpr):
    true_scores = []
    false_scores = []
    peaks = read_peaks(peaks_path)
    shuffled_peaks = creat_background(peaks, length, counter)
    write_fasta(shuffled_peaks, tmp_dir, "shuffled")
    write_fasta(peaks, tmp_dir, "train")
    make_inmode(peaks_path, path_to_inmode, path_to_java, length, order, tmp_dir, 'current')
    for true_score in true_scores_inmode(path_to_inmode, path_to_java, length, tmp_dir, "train", 'current'):
        true_scores.append(true_score)
    for false_score in false_scores_inmode(path_to_inmode, path_to_java, length, tmp_dir, "shuffled", 'current'):
        false_scores.append(false_score)
    roc = calculate_roc(true_scores, false_scores)
    fpr_current = fpr_at_tpr(true_scores, false_scores, tpr)
    auc_current = calculate_particial_auc(roc[0], roc[1], pfpr)
    print("Length {0}, Order {1};".format(length, order),
          "pAUC at {0} = {1};".format(pfpr, auc_current),
          "FPR = {0} at TPR = {1}".format(fpr_current, tpr))
    for length in range(length + 2, 34, 2):
        true_scores = []
        false_scores = []
        shuffled_peaks = creat_background(peaks, length, counter)
        write_fasta(shuffled_peaks, tmp_dir, "shuffled")
        make_inmode(peaks_path, path_to_inmode, path_to_java, length, order, tmp_dir, 'new')
        for true_score in true_scores_inmode(path_to_inmode, path_to_java, length, tmp_dir, "train", 'new'):
            true_scores.append(true_score)
        for false_score in false_scores_inmode(path_to_inmode, path_to_java, length, tmp_dir, "shuffled", 'new'):
            false_scores.append(false_score)
        roc = calculate_roc(true_scores, false_scores)
        fpr_new = fpr_at_tpr(true_scores, false_scores, tpr)
        auc_new = calculate_particial_auc(roc[0], roc[1], pfpr)
        print("Length {0}, Order {1};".format(length, order),
              "pAUC at {0} = {1};".format(pfpr, auc_new),
              "FPR = {0} at TPR = {1}".format(fpr_new, tpr))
        if auc_new > auc_current:
            shutil.copy(tmp_dir + '/new_inmode_model.xml',
                       tmp_dir + '/current_inmode_model.xml')
            auc_current = auc_new
        elif auc_new < auc_current and order <= 2:
            order += 1
        else:
            break
    return(order)


def de_novo_with_oprimization_inmode(peaks_path, length, path_to_inmode, 
    path_to_java, tmp_dir, output_path, tpr, pfpr):
    counter = 6000000
    order = 2
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    inmode_order = learn_optimized_inmode(peaks_path, counter, order, length, 
        path_to_inmode, path_to_java, tmp_dir, tpr, pfpr)
    shutil.copy(tmp_dir + '/new_inmode_model.xml', output_path)
    shutil.rmtree(tmp_dir)
    return(inmode_order)
