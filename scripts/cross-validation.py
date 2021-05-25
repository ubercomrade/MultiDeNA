import random
import math
import shutil
import os
import subprocess
import argparse
import sys
import functools
from multiprocessing import Pool
import numpy as np
from strum import strum
from operator import itemgetter
from lib.common import read_bamm, read_peaks, read_pwm, read_dipwm, \
calculate_merged_roc, calculate_short_roc, write_roc, \
calculate_fprs, creat_background, read_strum, calculate_particial_auc
from lib.speedup import score_dipwm, score_pwm, score_bamm


def complement(seq):
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


def write_fasta(sites, tmp_dir, tag):
    with open('{0}/{1}.fa'.format(tmp_dir, tag), 'w') as file:
        for index, site in enumerate(sites):
            file.write('>{0}\n{1}\n'.format(index, site))
    return(0)


# def get_motif_length(path):
#     with open(path, 'r') as file: ####
#         for i in file:
#             if i.startswith('>'):
#                 continue
#             else:
#                 motif_length = len(i.strip())
#                 break
#     file.close()
#     return(motif_length)


def get_model_length(auc_path):
    auc = []
    with open(auc_path) as file:
        for line in file:
            auc.append(tuple(map(float, line.strip().split())))
        file.close()
    auc.sort(key=itemgetter(1))
    best_auc = auc[-1][1]
    best_length = int(auc[-1][0])
    return(best_length)


def get_model_length_mm(auc_path):
    auc = []
    with open(auc_path) as file:
        for line in file:
            auc.append(tuple(map(float, line.strip().split())))
        file.close()
    auc.sort(key=itemgetter(2))
    best_auc = auc[-1][2]
    best_length = int(auc[-1][1])
    best_order = int(auc[-1][0])
    return(best_length)



def write_auc(path, auc):
    with open(path, 'w') as file:
        file.write('{0}'.format(auc))
    pass


def read_roc(path):
    container = dict()
    with open(path) as file:
        head = file.readline().strip().split()
        container[head[0]] = []
        container[head[1]] = []
        for line in file:
            line = line.strip().split()
            container[head[0]].append(float(line[0]))
            container[head[1]].append(float(line[1]))
    return(container)


def copy_results_of_cv(outdir, models_dir, model, pfpr):
    roc = read_roc(models_dir + '/{}_model/bootstrap.txt'.format(model))
    auc = calculate_particial_auc(roc['TPR'], roc['FPR'], pfpr)
    write_auc(outdir + '/{}_auc.txt'.format(model), auc)
    write_roc(outdir + '/{}_cv.txt'.format(model), roc)
    return(0)



#### PWM ####
def false_scores_pwm(peaks, pwm, length_of_site):
    false_scores = []
    append = false_scores.append
    for peak in peaks:
        complement_peak = complement(peak)
        full_peak = peak + 'N' * length_of_site + complement_peak
        n = len(full_peak) - length_of_site + 1
        for i in range(n):
            site = peak[i:length_of_site + i]
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


def cross_validation_pwm(pwm, length, peaks_path, counter, output_dir, pfpr):
    true_scores = []
    false_scores = []
    peaks = read_peaks(peaks_path)
    shuffled_peaks = creat_background(peaks, length, counter)
    for true_score in true_scores_pwm(peaks, pwm, length):
        true_scores.append(true_score)
    for false_score in false_scores_pwm(shuffled_peaks, pwm, length):
        false_scores.append(false_score)
    fprs = calculate_fprs(true_scores, false_scores)
    roc = calculate_short_roc(fprs, step=1)
    auc = calculate_particial_auc(roc['TPR'], roc['FPR'], pfpr)
    write_auc(output_dir + '/pwm_auc.txt', auc)
    write_roc(output_dir + "/pwm_cv.txt", roc)
    return(0)


#### INMODE ####
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


def true_scores_inmode(model_path, path_to_inmode, path_to_java, tmp_dir, fasta_tag, model_tag):
    scores = []
    args = [path_to_java, '-Xmx16G', '-Xms1G', 
            '-jar',
            path_to_inmode, 'scan',
            'i={0}'.format(model_path),
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


def false_scores_inmode(model_path, path_to_inmode, path_to_java, tmp_dir, fasta_tag, model_tag):
    scores = []
    args = [path_to_java, '-Xmx16G', '-Xms1G', 
            '-jar',
            path_to_inmode, 'scan',
            'i={0}'.format(model_path),
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


def cross_validation_inmode(inmode_model, peaks_path, counter, length, 
    path_to_inmode, path_to_java, tmp_dir, output_dir, pfpr):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    true_scores = []
    false_scores = []
    peaks = read_peaks(peaks_path)
    shuffled_peaks = creat_background(peaks, length, counter)
    write_fasta(shuffled_peaks, tmp_dir, "shuffled")
    write_fasta(peaks, tmp_dir, "train")
    for true_score in true_scores_inmode(inmode_model, path_to_inmode, path_to_java, tmp_dir, "train", 'current'):
        true_scores.append(true_score)
    for false_score in false_scores_inmode(inmode_model, path_to_inmode, path_to_java, tmp_dir, "shuffled", 'current'):
        false_scores.append(false_score)
    fprs = calculate_fprs(true_scores, false_scores)
    roc = calculate_short_roc(fprs, step=1)
    auc = calculate_particial_auc(roc['TPR'], roc['FPR'], pfpr)
    write_auc(output_dir + '/inmode_auc.txt', auc)
    write_roc(output_dir + "/inmode_cv.txt", roc)
    shutil.rmtree(tmp_dir)
    return(0)


#### BAMM ####
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


def cross_validation_bamm(bamm_path, bg_path, peaks_path, counter, length, output_dir, pfpr):
    true_scores = []
    false_scores = []
    bamm, order = read_bamm(bamm_path, bg_path)
    peaks = read_peaks(peaks_path)
    shuffled_peaks = creat_background(peaks, length, counter)
    for true_score in true_scores_bamm(peaks, bamm, order, length):
        true_scores.append(true_score)
    for false_score in false_scores_bamm(shuffled_peaks, bamm, order, length):
        false_scores.append(false_score)
    fprs = calculate_fprs(true_scores, false_scores)
    roc = calculate_short_roc(fprs, step=1)
    auc = calculate_particial_auc(roc['TPR'], roc['FPR'], pfpr)
    write_auc(output_dir + '/bamm_auc.txt', auc)
    write_roc(output_dir + "/bamm_cv.txt", roc)
    return(0)


#### DIPWM ####
def false_scores_dipwm(peaks, dipwm, length_of_site):
    false_scores = []
    append = false_scores.append
    for peak in peaks:
        complement_peak = complement(peak)
        full_peak = peak + 'N' * length_of_site + complement_peak
        n = len(full_peak) - length_of_site + 1
        for i in range(n):
            site = peak[i:length_of_site + i]
            if 'N' in site:
                continue
            score = score_dipwm(site, dipwm)
            false_scores.append(score)
    return(false_scores)


def true_scores_dipwm(peaks, dipwm, length_of_site):
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


def cross_validation_dipwm(dipwm, length, peaks_path, counter, output_dir, pfpr):
    true_scores = []
    false_scores = []
    peaks = read_peaks(peaks_path)
    shuffled_peaks = creat_background(peaks, length, counter)
    for true_score in true_scores_dipwm(peaks, dipwm, length):
        true_scores.append(true_score)
    for false_score in false_scores_dipwm(shuffled_peaks, dipwm, length):
        false_scores.append(false_score)
    fprs = calculate_fprs(true_scores, false_scores)
    roc = calculate_short_roc(fprs, step=1)
    auc = calculate_particial_auc(roc['TPR'], roc['FPR'], pfpr)
    write_auc(output_dir + '/dipwm_auc.txt', auc)
    write_roc(output_dir + "/dipwm_cv.txt", roc)
    return(0)


#### STRUM ####
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


def cross_validation_strum(strum_path, length, peaks_path, counter, output_dir, pfpr):
    true_scores = []
    false_scores = []
    peaks = read_peaks(peaks_path)
    strum_model = read_strum(strum_path)
    shuffled_peaks = creat_background(peaks, length, counter)
    for true_score, site in zip(*true_scores_strum(peaks, strum_model, length)):
        true_scores.append(true_score)
    for false_score in false_scores_strum(shuffled_peaks, strum_model):
        false_scores.append(false_score)
    fprs = calculate_fprs(true_scores, false_scores)
    roc = calculate_short_roc(fprs, step=1)
    auc = calculate_particial_auc(roc['TPR'], roc['FPR'], pfpr)
    write_auc(output_dir + '/strum_auc.txt', auc)
    write_roc(output_dir + "/strum_cv.txt", roc)
    return(0)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('data', action='store', help='directory contained MultiDeNa results')
    parser.add_argument('outdir', action='store', help='directory to write results')
    parser.add_argument('-i', '--ignore', action='store', help='names (short) of dir in <data> to ignore', type=list, required=False,
                       dest="ignore", default=[], nargs='+')
    parser.add_argument('-t', '--tools', action='store',  choices=['pwm', 'dipwm', 'bamm', 'inmode', 'sitega', 'strum'], metavar='N', nargs='+',
         help='list of models to use (pwm, dipwm, bamm, inmode, sitega, strum)', required=False,
         default=['pwm', 'bamm', 'inmode'])
    parser.add_argument('-I', '--inmode', action='store', dest='inmode',
                        required=True, help='path to inmode')
    parser.add_argument('-J', '--java', action='store', dest='java',
                    required=False, default="java", help='path to Java')
    parser.add_argument('-c', '--chipmunk', action='store', dest='chipmunk',
                        required=True, help='path to chipmunk')
    parser.add_argument('-C', '--processes', action='store', type=int, dest='cpu_count',
                        required=False, default=2, help='Number of processes to use, default: 2')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def run_cross_validation_through_all_data(t1, data_dir, write_dir, tags, tools,
    path_to_java, path_to_chipmunk, path_to_inmode):
    counter = 1000000
    pfpr = 0.001
    models = '{0}/{1}/models'.format(data_dir, t1)
    auc_dir = '{0}/{1}/auc'.format(data_dir, t1)
    pwm_model = read_pwm(models + '/pwm_model/pwm_model.pwm')
    dipwm_model = read_dipwm(models + '/dipwm_model/dipwm_model.pwm')
    inmode_model = models + '/inmode_model/inmode_model.xml'
    bamm_model = models + '/bamm_model/bamm_model.ihbcp'
    bg_bamm_model = models + '/bamm_model/bamm.hbcp'
    strum_model = models + '/strum_model/strum_model.pickle'
    for t2 in tags:
        fasta = '{0}/{1}/fasta'.format(data_dir, t2)
        peaks_path = fasta + '/train_sample.fa'
        results = '{0}/{1}/cv_{2}'.format(write_dir, t1, t2)
        if not os.path.isdir('{0}/{1}/'.format(write_dir, t1)):
            os.mkdir('{0}/{1}/'.format(write_dir, t1))
        if not os.path.isdir(results):
            os.mkdir(results)
        if 'pwm' in tools:
            if t1 != t2:
                print('{0} PWM model on data {1}'.format(t1, t2))
                motif_length = get_model_length(auc_dir + '/pwm/auc.txt')
                cross_validation_pwm(pwm_model, motif_length, peaks_path, counter, results, pfpr)
            else:
                copy_results_of_cv(results, models, 'pwm', pfpr)
                
        if 'dipwm' in tools:
            if t1 != t2:
                print('{0} diPWM model on data {1}'.format(t1, t2))
                motif_length = get_model_length(auc_dir + '/dipwm/auc.txt')
                cross_validation_dipwm(dipwm_model,  motif_length, peaks_path, counter, results, pfpr)
            else:
                copy_results_of_cv(results, models, 'dipwm', pfpr)

        if 'inmode' in tools:
            if t1 != t2:
                print('{0} InMoDe model on data {1}'.format(t1, t2))
                motif_length = get_model_length_mm(auc_dir + '/inmode/auc.txt')
                cross_validation_inmode(inmode_model, peaks_path, counter, motif_length, 
                                        path_to_inmode, path_to_java, results + '/inmode.tmp', results, pfpr)
            else:
                copy_results_of_cv(results, models, 'inmode', pfpr)

        if 'bamm' in tools:
            if t1 != t2:
                print('{0} BaMM model on data {1}'.format(t1, t2))
                motif_length = get_model_length_mm(auc_dir + '/bamm/auc.txt')
                cross_validation_bamm(bamm_model, bg_bamm_model, peaks_path, counter, motif_length, results, pfpr)
            else:
                copy_results_of_cv(results, models, 'bamm', pfpr)

        if 'strum' in tools:
            if t1 != t2:
                print('{0} StruM model on data {1}'.format(t1, t2))
                motif_length = get_model_length(auc_dir + '/strum/auc.txt')
                cross_validation_strum(strum_model, motif_length, peaks_path, counter, results, pfpr)
            else:
                copy_results_of_cv(results, models, 'strum', pfpr)
    pass


def main():
    # wd = "/Users/tsukanov/Documents/PhD/remap/foxa2/results/"
    # write_dir = "//Users/tsukanov/Documents/PhD/remap/foxa2/cross-validation/"
    # ignore = ['ENCSR066EBK.FOXA2.Hep-G2', 'GSE90454.FOXA2.BJ1-hTERT_FoxHnf1aCoExp', 'GSE90454.FOXA2.BJ1-hTERT_GATA4',
    #           'GSE90454.FOXA2.BJ1-hTERT_Mimo', 'GSE90454.FOXA2.BJ1-hTERT_Mimo_Release', 'GSE90454.FOXA2.Hep-G2']
    args = parse_args()
    data_dir = args.data
    write_dir = args.outdir
    ignore = args.ignore
    tools = args.tools

    path_to_java = args.java
    path_to_chipmunk = args.chipmunk
    path_to_inmode = args.inmode
    cpu_count = args.cpu_count
    
    tags = [i for i in os.listdir(data_dir) if not i.startswith('.') and not i in ignore]
    # for t1 in tags:
    #     run_cross_validation_through_all_data(t1, data_dir, write_dir, tags, tools,
    #         path_to_java, path_to_chipmunk, path_to_inmode)
    with Pool(cpu_count) as p:
        p.map(functools.partial(run_cross_validation_through_all_data, data_dir=data_dir, 
            write_dir=write_dir, 
            tags=tags, 
            tools=tools,
            path_to_java=path_to_java, 
            path_to_chipmunk=path_to_chipmunk, 
            path_to_inmode=path_to_inmode), tags)
    return(0)

if __name__=="__main__":
    main()