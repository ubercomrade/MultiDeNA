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
from multidena.lib.common import calculate_particial_auc, \
write_auc_with_order, calculate_merged_roc, calculate_short_roc, \
write_roc, calculate_fprs, write_table


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


def best_scores_inmode(path_to_inmode, path_to_java, motif_length, tmp_dir, fasta_tag, model_tag, order):
    scores = []
    args = [path_to_java, '-Xmx16G', '-Xms1G', 
            '-jar',
            path_to_inmode, 'scan',
            'i={0}/inmode_model_{1}_{2}.xml'.format(tmp_dir, order, model_tag),
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


def false_scores_inmode(path_to_inmode, path_to_java, motif_length, tmp_dir, fasta_tag, model_tag, order):
    scores = []
    args = [path_to_java, '-Xmx16G', '-Xms1G', 
            '-jar',
            path_to_inmode, 'scan',
            'i={0}/inmode_model_{1}_{2}.xml'.format(tmp_dir, order, model_tag),
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
                tmp_dir + '/inmode_model_{0}_{1}.xml'.format(order, tag))
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


def learn_optimized_inmode(peaks_path, backgroud_path, counter, 
    path_to_inmode, path_to_java, tmp_dir, output_auc, pfpr):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    if not os.path.isdir(output_auc):
        os.mkdir(output_auc)
    open(output_auc + '/auc.txt', 'w').close()
    for order in range(1,5):
        #for length in range(12, 41, 4):
        for length in range(8, 21, 4):
            true_scores = []
            false_scores_roc = []
            false_scores_prc = []
            peaks = read_peaks(peaks_path)
            for step in ['odd', 'even']:
                if step == 'odd':
                    train_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 != 0]
                    test_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 == 0]
                else:
                    train_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 == 0]
                    test_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 != 0]                
                if os.path.isfile(backgroud_path):
                    shuffled_peaks = read_peaks(backgroud_path)
                else:
                    shuffled_peaks = creat_background(test_peaks, length, counter)
                write_fasta(shuffled_peaks, tmp_dir, "shuffled")
                write_fasta(train_peaks, tmp_dir, "train")
                write_fasta(test_peaks, tmp_dir, "test")
                make_inmode('{0}/{1}.fa'.format(tmp_dir, 'train'), path_to_inmode, path_to_java, length, order, tmp_dir, str(length))
                for true_score in best_scores_inmode(path_to_inmode, path_to_java, length, tmp_dir, "test", str(length), str(order)):
                    true_scores.append(true_score)
                for false_score in best_scores_inmode(path_to_inmode, path_to_java, length, tmp_dir, "shuffled", str(length), str(order)):
                    false_scores_prc.append(true_score)
                for false_score in false_scores_inmode(path_to_inmode, path_to_java, length, tmp_dir, "shuffled", str(length), str(order)):
                    false_scores_roc.append(false_score)
                shutil.copy(tmp_dir + '/inmode_model_{0}_{1}.xml'.format(order, length),
                            output_auc + '/inmode_model_{0}_{1}_{2}.xml'.format(step, order, length))            
            fprs = calculate_fprs(true_scores, false_scores_roc)
            prc = calculate_prc(true_scores, false_scores_prc)
            roc = calculate_short_roc(fprs, step=1)
            merged_roc = calculate_merged_roc(fprs)
            auc_roc = calculate_particial_auc(merged_roc['TPR'], merged_roc['FPR'], pfpr)
            auc_prc = calculate_particial_auc(prc['PRECISION'], prc['RECALL'], 1.01)
            print("Length {0}; Order {1}".format(length, order), "pAUC at {0} = {1};".format(pfpr, auc_roc), "PRC = {0}".format(auc_prc))
            write_table(output_auc + '/statistics.txt', auc_roc, auc_prc, order, length)
            write_roc(output_auc + "/training_roc_merged_{0}_{1}.txt".format(length, order), merged_roc)
            write_roc(output_auc + "/training_roc_{0}_{1}.txt".format(length, order), roc)
            write_roc(output_auc + "/training_prc_{0}_{1}.txt".format(length, order), prc) 
    shutil.rmtree(tmp_dir)
    return(0)


def choose_best_model(output_auc):
    auc = []
    with open(output_auc + '/statistics.txt') as file:
        for line in file:
            auc.append(tuple(map(float, line.strip().split())))
        file.close()
    auc.sort(key=itemgetter(-1))
    order = int(auc[-1][1])
    length = int(auc[-1][0])
    return(length, order)


def de_novo_with_oprimization_inmode(peaks_path, backgroud_path, path_to_inmode, 
    path_to_java, tmp_dir, output_dir, output_auc, pfpr):
    counter = 1000000
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    else:
        shutil.rmtree(tmp_dir)
        os.mkdir(tmp_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.exists(output_auc):
        os.mkdir(output_auc)
    if os.path.exists(output_auc + '/auc.txt'):
        os.remove(output_auc + '/auc.txt')
    learn_optimized_inmode(peaks_path, backgroud_path, counter,
                           path_to_inmode, path_to_java, 
                           tmp_dir, output_auc, pfpr)
    length, order = choose_best_model(output_auc)
    
    shutil.copy(output_auc + '/training_prc_{0}_{1}.txt'.format(length, order), 
             output_dir + '/prc.txt')
    shutil.copy(output_auc + '/training_roc_{0}_{1}.txt'.format(length, order), 
             output_dir + '/roc.txt')
    shutil.copy(output_auc + '/training_roc_merged_{0}_{1}.txt'.format(length, order), 
             output_dir + '/roc_merged.txt')
    
    args = [path_to_java, '-Xmx16G', '-Xms1G', '-jar', path_to_inmode,
    'denovo', 'i={}'.format(peaks_path), 'm={}'.format(length), 'outdir={}'.format(tmp_dir),
    'mo={}'.format(order)]
    r = subprocess.run(args, capture_output=True)
    shutil.copy(tmp_dir + '/Learned_DeNovo({0},{1},2)_motif/XML_of_DeNovo({0},{1},2)_motif.xml'.format(length,order),
                output_dir + '/inmode_model.xml')
    shutil.rmtree(tmp_dir)
    with open(f"{output_dir}/properties.txt", "w") as file:
        file.write(f"{length}\t{order}\n")
    return 0
