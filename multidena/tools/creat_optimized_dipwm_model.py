import random
import shutil
import os
import subprocess
import bisect
from operator import itemgetter
from shutil import copyfile
from multidena.lib.common import read_peaks, sites_to_dipwm, creat_background, \
write_fasta, complement, make_dipcm, make_dipfm, \
make_dipwm, write_dipwm, write_dipfm, \
calculate_particial_auc, write_auc, \
calculate_merged_roc, calculate_short_roc, \
write_roc, calculate_fprs, write_table, calculate_prc
from multidena.lib.speedup import creat_table_bootstrap, score_dipwm


def run_di_chipmunk(path_to_java, path_to_chipmunk, fasta_path, path_out, motif_length_start, motif_length_end, cpu_count):
    args = [path_to_java, '-cp', path_to_chipmunk,
                   'ru.autosome.di.ChIPMunk', str(motif_length_start), str(motif_length_end), 'yes', '1.0',
                   's:{}'.format(fasta_path),
                  '100', '10', '1', str(cpu_count), 'random']
    p = subprocess.run(args, shell=False, capture_output=True)
    out = p.stdout
    with open(path_out, 'wb') as file:
        file.write(out)
    return(0)


def parse_chipmunk(path):
    with open(path, 'r') as file:
        container = []
        for line in file:
            d = {'name': str(), 'start': int(), 'end': int(),
                 'seq': str(), 'strand': str()}
            if line.startswith('WORD|'):
                line = line[5:].strip()
                line = line.split()
                d['name'] = 'peaks_' + str(int(line[0]) - 1)
                d['start'] = int(line[1])
                d['end'] = int(line[1]) + len(line[2])
                d['seq'] = line[2]
                if line[4] == 'direct':
                    d['strand'] = '+'
                else:
                    d['strand'] = '-'
                container.append(d)
            else:
                continue
    seqs = [i['seq'] for i in container if len(set(i['seq']) - {'A', 'C', 'G', 'T'}) == 0]
    return(seqs)


def false_scores_dipwm(peaks, dipwm, length_of_site):
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
            score = score_dipwm(site, dipwm)
            false_scores.append(score)
    return(false_scores)


def best_scores_dipwm(peaks, dipwm, length_of_site):
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
            score = score_dipwm(site, dipwm)
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


def learn_optimized_dipwm(peaks_path, backgroud_path, counter, path_to_java, path_to_chipmunk, tmp_dir, output_auc, cpu_count, pfpr, min_length, max_length, length_step):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    if not os.path.isdir(output_auc):
        os.mkdir(output_auc)
    if os.path.exists(output_auc + '/auc.txt'):
        os.remove(output_auc + '/auc.txt')
    #for length in range(12, 41, 4):
    for length in range(min_length, max_length + 1, length_step):
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
            write_fasta(train_peaks, tmp_dir + '/train.fasta')
            if os.path.isfile(backgroud_path):
                shuffled_peaks = read_peaks(backgroud_path)
            else:
                shuffled_peaks = creat_background(test_peaks, length, counter)
            run_di_chipmunk(path_to_java, path_to_chipmunk,
                         tmp_dir + '/train.fasta', tmp_dir + '/chipmunk_results.txt',
                         length, length, cpu_count)
            sites = parse_chipmunk(tmp_dir + '/chipmunk_results.txt')
            sites = list(set(sites))
            dipwm = sites_to_dipwm(sites)
            for true_score in best_scores_dipwm(peaks, dipwm, length):
                true_scores.append(true_score)
            for false_score in best_scores_dipwm(shuffled_peaks, dipwm, length):
                false_scores_prc.append(false_score)
            for false_score in false_scores_dipwm(shuffled_peaks, dipwm, length):
                false_scores_roc.append(false_score)
            dipcm = make_dipcm(sites)
            dipfm = make_dipfm(dipcm)
            dipwm = make_dipwm(dipfm)
            tag = 'dipwm_model_{0}_{1}'.format(step, length)
            write_dipwm(output_auc, tag, dipwm)
            write_dipfm(output_auc, tag, dipfm)
            write_sites(output=output_auc, tag=tag, sites=sites)
        fprs = calculate_fprs(true_scores, false_scores_roc)
        prc = calculate_prc(true_scores, false_scores_prc)
        roc = calculate_short_roc(fprs, step=1)
        merged_roc = calculate_merged_roc(fprs)
        auc_roc = calculate_particial_auc(merged_roc['TPR'], merged_roc['FPR'], pfpr)
        auc_prc = calculate_particial_auc(prc['PRECISION'], prc['RECALL'], 1.01)
        print("Length {};".format(length), "pAUC at {0} = {1};".format(pfpr, auc_roc), "PRC = {0}".format(auc_prc))
        write_table(output_auc + '/statistics.txt', auc_roc, auc_prc, length)
        write_roc(output_auc + "/training_roc_merged_{0}.txt".format(length), merged_roc)
        write_roc(output_auc + "/training_roc_{0}.txt".format(length), roc)
        write_roc(output_auc + "/training_prc_{0}.txt".format(length), prc)
    shutil.rmtree(tmp_dir)
    return(0)


def choose_best_model(output_auc):
    auc = []
    with open(output_auc + '/statistics.txt') as file:
        for line in file:
            auc.append(tuple(map(float, line.strip().split())))
        file.close()
    auc.sort(key=itemgetter(-1))
    length = int(auc[-1][0])
    return(length)


def de_novo_with_oprimization_dipwm(peaks_path, backgroud_path, path_to_java, path_to_chipmunk,
    tmp_dir, output_dir, output_auc, cpu_count, pfpr, min_length, max_length, length_step):
    counter = 1000000
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    learn_optimized_dipwm(peaks_path, backgroud_path, counter, path_to_java,
       path_to_chipmunk, tmp_dir, output_auc, cpu_count, pfpr, min_length, max_length, length_step)
    length = choose_best_model(output_auc)
    copyfile(output_auc + '/training_prc_{}.txt'.format(length),
             output_dir + '/prc.txt')
    copyfile(output_auc + '/training_roc_{}.txt'.format(length),
             output_dir + '/roc.txt')
    copyfile(output_auc + '/training_roc_merged_{}.txt'.format(length),
             output_dir + '/roc_merged.txt')
    run_di_chipmunk(
        path_to_java, path_to_chipmunk,
        peaks_path,
        output_dir + '/chipmunk_results.txt',
        length, length, cpu_count)
    sites = parse_chipmunk(output_dir + '/chipmunk_results.txt')
    dipcm = make_dipcm(sites)
    dipfm = make_dipfm(dipcm)
    dipwm = make_dipwm(dipfm)
    nsites = len(sites)
    background = {'A': 0.25,
                 'C': 0.25,
                 'G': 0.25,
                 'T': 0.25}
    tag = 'dipwm_model'
    write_dipwm(output_dir, tag, dipwm)
    write_dipfm(output_dir, tag, dipfm)
    write_sites(output=output_dir, tag=tag, sites=sites)
    return(length)
