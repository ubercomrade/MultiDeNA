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


def read_peaks(path):
    container = []
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                append(line.strip().upper())
    return(container)


def sitega_bootstrap(data_dir, lpd_length=6, start_motif_length=8, end_motif_length=20, step=4):
    args = ['andy0bsn5',  f'{data_dir}/',  'train.fa',  'background.fa',
    f'{lpd_length}', f'{start_motif_length}', f'{end_motif_length}',
    f'{step}', '-1', '2', '6', f'{data_dir}/']
    print(args)
    capture = subprocess.run(args, capture_output=False)
    return 0


def get_sitega_model(output_dir, tmp_dir, motif_length, number_of_lpds):
    args = ['andy05', f'{tmp_dir}/', 
    'train.fa', 'background.fa', 
    '6', f'{motif_length}', 
    f'{number_of_lpds}', '6', f'{tmp_dir}/']
    print(args)
    capture = subprocess.run(args, capture_output=False)
    shutil.copy(f'{tmp_dir}/train.fa_mat', f'{output_dir}/sitega.mat')
    return 0
    
    
def parse_roc(path_in, path_out):
    container = []
    with open(path_in) as file:
        container.append(file.readline().strip().split())
        container.append(file.readline().strip().split()[1:])
        file.close()
    with open(path_out, 'w') as file:
        file.write('TPR\tFPR\n')
        for line in zip(*container):
            file.write('\t'.join(line) + '\n')
    return 0


def parse_auc_table(path_in, path_out):
    container = []
    with open(path_in) as file:
        for line in file:
            line = line.strip().split()[1:]
            container.append(line)
        file.close()
    with open(path_out, 'w') as file:
        for line in container:
            file.write(f'{line[1]}\t{line[0]}\t{line[2]}\n')
        file.close()
    return 0

    
def write_fasta(sites, tmp_dir, tag):
    with open('{0}/{1}.fa'.format(tmp_dir, tag), 'w') as file:
        for index, site in enumerate(sites):
            file.write('>{0}\n{1}\n'.format(index, site))
    return(0)
    
    
def creat_background(peaks, counter):
    shuffled_peaks = []
    number_of_sites = 0
    length_of_site = 20
    while counter > number_of_sites:
        peak = random.choice(peaks)
        shuffled_peak = ''.join(random.sample(peak, len(peak)))
        shuffled_peaks.append(shuffled_peak)
        number_of_sites += (len(''.join(shuffled_peak)) - length_of_site + 1) * 2
    return(shuffled_peaks)


def complement(seq):
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


def learn_optimized_sitega(peaks_path, backgroud_path, counter, tmp_dir, output_dir, output_auc):
    if os.path.isfile(backgroud_path):
        shutil.copy(backgroud_path, f'{tmp_dir}/background.fa')
    else:
        shuffled_peaks = creat_background(peaks_path, counter)
        write_fasta(shuffled_peaks, tmp_dir, 'background')
    shutil.copy(peaks_path, f'{tmp_dir}/train.fa')
    sitega_bootstrap(tmp_dir, lpd_length=6, start_motif_length=8, end_motif_length=20, step=4)
    shutil.copy(f'{tmp_dir}/train.fa_mat', f'{output_auc}/sitega_bootstrap_models.fa_mat')
    shutil.copy(f'{tmp_dir}/train.fa_roc_bs.txt', f'{output_auc}/sitega_bootstrap_rocs.txt')
    parse_auc_table(f'{tmp_dir}/train.fa_auc_bs.txt',
              f'{output_auc}/auc.txt')
    parse_roc(f'{tmp_dir}/train.fa_best_roc_bs.txt',
              f'{output_dir}/bootstrap.txt')    
    return(0)


def choose_best_model(output_auc):
    auc = []
    with open(output_auc + '/auc.txt') as file:
        for line in file:
            auc.append(tuple(map(float, line.strip().split())))
        file.close()
    auc.sort(key=itemgetter(-1))
    order = int(auc[-1][0])
    length = int(auc[-1][1])
    return(length, order)


def de_novo_with_oprimization_sitega(peaks_path, backgroud_path, 
                                     tmp_dir, output_dir, output_auc, pfpr=0.001):
    counter = 5000000
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
    learn_optimized_sitega(peaks_path, backgroud_path, counter, 
                           tmp_dir, output_dir, output_auc)
    motif_length, number_of_lpds = choose_best_model(output_auc)
    get_sitega_model(output_dir, tmp_dir, motif_length, number_of_lpds)
    shutil.rmtree(tmp_dir)
    return motif_length, number_of_lpds

