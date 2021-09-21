import random
import shutil
import os
import subprocess
import bisect
from shutil import copyfile
from operator import itemgetter
from lib.common import read_peaks, sites_to_pwm, creat_background, \
write_fasta, complement, make_pcm, make_pfm, \
make_pwm, write_pwm, write_pfm, write_meme, \
calculate_particial_auc, write_auc, calculate_merged_roc, \
calculate_short_roc, write_roc, calculate_fprs
from lib.speedup import creat_table_bootstrap, score_pwm


def run_streme(fasta_path, backgroud_path, dir_out, motif_length):
    args = ['streme', '--p', fasta_path,
            '--n', backgroud_path,
           '--oc', dir_out,
           '--objfun', 'de',
           '--w', str(motif_length),
           '-nmotifs', '3']
    #print(' '.join(args))
    p = subprocess.run(args, shell=False, capture_output=True)
    #print(p.stdout)
    return(0)


def run_streme_hmm_background(fasta_path, dir_out, motif_length):
    args = ['streme', '--p', fasta_path,
            '--kmer', '4',
           '--oc', dir_out,
           '--objfun', 'de',
           '--w', str(motif_length),
           '-nmotifs', '3']
    #print(' '.join(args))
    p = subprocess.run(args, shell=False, capture_output=True)
    #print(p.stdout)
    return(0)


def parse_streme(path):
    pfm = {'A': [], 'C': [], 'G': [], 'T': []}
    with open(path) as file:
        for line in file:
            if line.startswith('Background'):
                line = file.readline().strip().split()
                background = {'A': float(line[1]),
                             'C': float(line[3]),
                             'G': float(line[5]),
                             'T': float(line[7])}
            elif line.startswith('letter-probability'):
                break
        line = line.strip().split()
        length = int(line[5])
        nsites = int(line[7])
        for i in range(length):
            line = file.readline().strip().split()
            for letter, value in zip(pfm.keys(), line):
                    pfm[letter].append(float(value))
    file.close()
    return pfm, background, length, nsites


def false_scores_pwm(peaks, pwm, length_of_site):
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
            score = score_pwm(site, pwm)
            false_scores.append(score)
    return(false_scores)


def true_scores_pwm(peaks, pwm, length_of_site):
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
            score = score_pwm(site, pwm)
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


def write_sites(output, tag, sites):
    with open(output + '/' + tag + '.fasta', 'w') as file:
        for index, site in enumerate(sites):
            file.write(site + '\n')
    return(0)


def learn_optimized_pwm(peaks_path, backgroud_path, counter, tmp_dir, output_auc, pfpr):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    if not os.path.isdir(output_auc):
        os.mkdir(output_auc)
    if os.path.exists(output_auc + '/auc.txt'):
        os.remove(output_auc + '/auc.txt')
    #for length in range(12, 41, 4):
    for length in range(10, 31, 4):
        true_scores = []
        false_scores = []
        peaks = read_peaks(peaks_path)
        for step in ['odd', 'even']:
            if step == 'odd':
                train_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 != 0]
                test_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 == 0]
            else:
                train_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 == 0]
                test_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 != 0]
            write_fasta(train_peaks, tmp_dir + '/train.fasta')
            if os.path.isfile(backgroud_path):
                shuffled_peaks = read_peaks(backgroud_path)
                run_streme(tmp_dir + '/train.fasta', backgroud_path, tmp_dir, length)
            else:
                shuffled_peaks = creat_background(test_peaks, length, counter)
                run_streme_hmm_background(tmp_dir + '/train.fasta', tmp_dir, length)
            pfm, background, length, nsites = parse_streme(tmp_dir + '/streme.txt')
            pwm = make_pwm(pfm)
            for true_score in true_scores_pwm(test_peaks, pwm, length):
                true_scores.append(true_score)
            for false_score in false_scores_pwm(shuffled_peaks, pwm, length):
                false_scores.append(false_score)
            tag = 'pwm_model_{0}_{1}'.format(step, length)
            write_meme(output_auc, tag, pfm, background, nsites)
            write_pwm(output_auc, tag, pwm)
            write_pfm(output_auc, tag, pfm)
        fprs = calculate_fprs(true_scores, false_scores)
        roc = calculate_short_roc(fprs, step=1)
        merged_roc = calculate_merged_roc(fprs)
        auc = calculate_particial_auc(merged_roc['TPR'], merged_roc['FPR'], pfpr)
        print("Length {};".format(length), "pAUC at {0} = {1};".format(pfpr, auc))
        write_auc(output_auc + '/auc.txt', auc, length)
        write_roc(output_auc + "/training_bootstrap_merged_{0}.txt".format(length), merged_roc)
        write_roc(output_auc + "/training_bootstrap_{0}.txt".format(length), roc)
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


def de_novo_with_oprimization_pwm_streme(peaks_path, backgroud_path, tmp_dir, output_dir, output_auc, pfpr):
    counter = 1000000
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    if not os.path.isdir(output_auc):
        os.mkdir(output_auc)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    learn_optimized_pwm(peaks_path, backgroud_path, counter, tmp_dir, output_auc, pfpr)
    length = choose_best_model(output_auc)
    copyfile(output_auc + '/training_bootstrap_{}.txt'.format(length), 
             output_dir + '/bootstrap.txt')
    copyfile(output_auc + '/training_bootstrap_merged_{}.txt'.format(length), 
             output_dir + '/bootstrap_merged.txt')
    run_streme(peaks_path, backgroud_path, output_dir, length)
    pfm, background, length, nsites = parse_streme(output_dir + '/streme.txt')
    pwm = make_pwm(pfm)
    tag = 'pwm_model'
    write_meme(output_dir, tag, pfm, background, nsites)
    write_pwm(output_dir, tag, pwm)
    write_pfm(output_dir, tag, pfm)
    return(length)
