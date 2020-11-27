import subprocess
import os
import random
import math
import shutil
import fnmatch

import functools
from multiprocessing import Pool
from operator import itemgetter
from tools.clear_from_n import clear_from_n
from lib.common import read_peaks, write_table_bootstrap, \
creat_background, complement, write_table_bootstrap_wide, \
calculate_roc
from lib.speedup import creat_table_bootstrap


def write_fasta(sites, tmp_dir, tag):
    with open('{0}/{1}.fa'.format(tmp_dir, tag), 'w') as file:
        for index, site in enumerate(sites):
            file.write('>{0}\n{1}\n'.format(index, site))
    return(0)


def best_scores_sitega(tmp_dir, tag):
    scores = []
    args = ['andy1_mat',
        '{}'.format(tmp_dir + '/{}.fa'.format(tag)),
        '{}'.format(tmp_dir + '/sitega.mat'),
        '{}'.format(tmp_dir + '/thr_table.txt'),
        '{}'.format(0),
       '{}'.format(tmp_dir + '/sitega_true.pro')]
    r = subprocess.run(args, capture_output=False)
    with open(tmp_dir + '/{}.fa_bestscosg'.format(tag)) as file:
        for line in file:
            scores.append(float(line.strip()))
    file.close()
    return(scores)


def make_sitega(tmp_dir, length, lpd):
    args = ['monte0dg' ,'6', tmp_dir + '/train.fa', tmp_dir + '/peaks.mnt']
    capture = subprocess.run(args, capture_output=True)
    args = ['andy02', tmp_dir + '/peaks.mnt', str(length), str(lpd), str(lpd), '10']
    capture = subprocess.run(args, capture_output=True)
    for file in os.listdir(tmp_dir):
        if fnmatch.fnmatch(file, 'train.fa_mat_*'):
            shutil.move('{0}/{1}'.format(tmp_dir, file), '{0}/sitega.mat'.format(tmp_dir))
    return(0)


def bootstrap_sitega(peaks, length_of_site, lpd, counter, tmp_dir):
    true_scores = []
    false_scores = []
    number_of_peaks = len(peaks)
    for i in range(5):
        print(i)
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)
            with open(tmp_dir + '/thr_table.txt', 'w') as file:
                file.write("0.0\t0.0")
            file.close()
        train_peaks = random.choices(peaks, k=round(0.9 * number_of_peaks))
        test_peaks = [peak for peak in peaks  if not peak in train_peaks]
        shuffled_peaks = creat_background(test_peaks, length_of_site, counter / 5)
        write_fasta(train_peaks, tmp_dir, "train")
        write_fasta(test_peaks, tmp_dir, "test")
        write_fasta(shuffled_peaks, tmp_dir, "shuffled")
        make_sitega(tmp_dir, length_of_site, lpd)
        for true_score in best_scores_sitega(tmp_dir, "test"):
            true_scores.append(true_score)
        print(true_scores[-10:])
        for false_score in best_scores_sitega(tmp_dir, "shuffled"):
            false_scores.append(false_score)
        print(false_scores[-10:])
        shutil.rmtree(tmp_dir)
    table = creat_table_bootstrap(true_scores, false_scores)
    table_full = calculate_roc(true_scores, false_scores)
    return(table, table_full)


# MULTIPROCESS WORK IN PROGRESS (NEDEED???)
def support(tmp_dir, peaks, length_of_site, lpd, counter):
    true_scores = []
    false_scores = []
    number_of_peaks = len(peaks)
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
        with open(tmp_dir + '/thr_table.txt', 'w') as file:
            file.write("0.0\t0.0")
        file.close()
    train_peaks = random.choices(peaks, k=round(0.9 * number_of_peaks))
    test_peaks = [peak for peak in peaks  if not peak in train_peaks]
    shuffled_peaks = creat_background(test_peaks, length_of_site, counter / 5)
    write_fasta(train_peaks, tmp_dir, "train")
    write_fasta(test_peaks, tmp_dir, "test")
    write_fasta(shuffled_peaks, tmp_dir, "shuffled")
    make_sitega(tmp_dir, length_of_site, lpd)
    for true_score in best_scores_sitega(tmp_dir, "test"):
        true_scores.append(true_score)
    for false_score in best_scores_sitega(tmp_dir, "shuffled"):
        false_scores.append(false_score)
    shutil.rmtree(tmp_dir)
    return(true_scores, false_scores)


def bootstrap_sitega_multiprocessing(peaks, length_of_site, lpd, counter, tmp_dir):
    true_scores = []
    false_scores = []
    if tmp_dir[-1] == '/':
        tmp_dir = tmp_dir[:-1]
    tmp_dirs = ["{0}_{1}".format(tmp_dir, i) for i in range(5)]
    with Pool(5) as p:
        results = p.map(functools.partial(support, peaks=peaks, length_of_site=length_of_site, lpd=lpd, counter=counter),
              tmp_dirs)
    true_scores = [true_score for attempt in results for true_score in attempt[0]]
    false_scores = [false_score for attempt in results for false_score in attempt[1]]
    table = creat_table_bootstrap(true_scores, false_scores)
    table_full = calculate_roc(true_scores, false_scores)
    return(table, table_full)


def bootstrap_for_sitega(peaks_path_no_n, results_path, results_path_wide, length_of_site, lpd, tmp_dir, counter=5000000):
    peaks = read_peaks(peaks_path_no_n)
    #table, table_full = bootstrap_sitega(peaks, length_of_site, lpd, counter, tmp_dir)
    table, table_full = bootstrap_sitega_multiprocessing(peaks, length_of_site, lpd, counter, tmp_dir)
    write_table_bootstrap(results_path, table)
    write_table_bootstrap_wide(results_path_wide, table_full)
    return(0)