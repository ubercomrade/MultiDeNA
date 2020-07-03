import math
import sys
import os
import random
import itertools
import argparse
from operator import itemgetter
from lib.common import make_pcm, make_pfm, make_pwm, sites_to_pwm, \
calculate_scores_pwm_bootstrap, write_meme, write_pwm, write_pfm


def parse_chipmunk_words(path):
    with open(path, 'r') as file:
        output = []
        for line in file:
            d = {'name': str(), 'start': int(), 'end': int(),
                 'seq': str(), 'strand': str()}
            if line.startswith('WORD|'):
                line = line[5:].strip()
                line = line.split()
                d['name'] = int(line[0])
                d['start'] = int(line[1])
                d['end'] = int(line[1]) + len(line[2])
                d['seq'] = line[2]
                if line[4] == 'direct':
                    d['strand'] = '+'
                else:
                    d['strand'] = '-'
                output.append(d)
            else:
                continue
    return(output)


def chipmunk_sites(fasta, chipmunk):
    motifs = []
    for line in chipmunk:
        if line['strand'] == '+':
            start = line['start']
            end = line['end']
            if start < 0:
                continue
            elif end >= len(fasta[line['name']]):
                continue
            else:
                motifs.append(fasta[line['name']][start:end])
        else:
            start = len(fasta[line['name']]) - line['end']
            end = len(fasta[line['name']]) - line['start']
            if start < 0:
                continue
            elif end >= len(fasta[line['name']]):
                continue
            else:
                motifs.append(complement(fasta[line['name']])[start:end])
    motifs = [i for i in motifs if not 'N' in i]
    return(motifs)


def read_fasta(path):
    sequences = []
    with open(path, 'r') as file:
        sequences = [i.strip().upper() for i in file if i.strip()[0] != '>']
    return(sequences)


def complement(seq):
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


def remove_equalent_seq(seq_list, homology=0.95):
    seq_list = list(seq_list)
    treshold = homology * len(seq_list[0])
    for seq1 in tuple(seq_list):
        sub_seq_list = list(seq_list)
        sub_seq_list.remove(seq1)
        for seq2 in sub_seq_list:
            score = len([i for i, j in zip(seq1, seq2) if i == j])
            if score >= treshold:
                seq_list.remove(seq1)
                break
    return(seq_list)


def calculate_fpr(sites, size_of_random_sample):
    true_scores = []
    false_scores = []
    tpr = 0.5
    number_of_sites = len(sites)
    len_of_site = len(sites[0])
    for i in range(10):
        train_sample = random.choices(sites, k=round(0.9 * number_of_sites))
        test_sample = [site for site in sites  if not site in train_sample]
        random_sample = [''.join(random.sample(list(random.choice(test_sample)), len_of_site)) for i in range(len(test_sample) * size_of_random_sample)]
        pwm = sites_to_pwm(train_sample)
        for true_score in calculate_scores_pwm_bootstrap(test_sample, pwm):
            true_scores.append(true_score)
        for false_score in calculate_scores_pwm_bootstrap(random_sample, pwm):
            false_scores.append(false_score)

    true_scores.sort(reverse=True)
    false_scores.sort(reverse=True)
    false_length = len(false_scores)
    true_length = len(true_scores)
    score = true_scores[round(true_length * tpr) - 1]
    fpr = sum([1 if false_score >= score else 0 for false_score in false_scores]) / false_length
    return(fpr)


def write_sites(output, tag, sites):
    with open(output + '/' + tag + '.fasta', 'w') as file:
        for index, site in enumerate(sites):
            file.write('>site_' + str(index) + '\n')
            file.write(site + '\n')


def make_optimized_pwm(chipmunk_path, fasta_path, output_dir, times, tag, cpu_count):
    chipmunk = parse_chipmunk_words(chipmunk_path)
    fasta = read_fasta(fasta_path)
    background = {'A': 0.25,
                 'C': 0.25,
                 'G': 0.25,
                 'T': 0.25}
        fprs = []
    optimal_sites = chipmunk_sites(fasta, chipmunk)
    #optimal_sites = remove_equalent_seq(optimal_sites, homology=0.95)
    optimal_fpr = calculate_fpr(optimal_sites, 500)
    fprs.append(optimal_fpr)
    print(optimal_fpr)

    for i in chipmunk:
        i['start'] -= 1
        i['end'] += 1
    sites = chipmunk_sites(fasta, chipmunk)
    #sites = remove_equalent_seq(sites, homology=0.95)
    fpr = calculate_fpr(sites, 500)
    fprs.append(fpr)
    print(fpr)

    while fpr <= optimal_fpr:
        optimal_fpr = fpr
        optimal_sites = sites.copy()
        for i in chipmunk:
            i['start'] -= 1
            i['end'] += 1
        sites = chipmunk_sites(fasta, chipmunk)
        #sites = remove_equalent_seq(sites, homology=0.95)
        fpr = calculate_fpr(sites, 500)
        fprs.append(fpr)
        if fpr <= optimal_fpr and 1 - fprs[-1]/fprs[-2] < 0.1:
            print(fpr)
            break
        print(fpr)

    pcm = make_pcm(sites)
    pfm = make_pfm(pcm)
    pwm = make_pwm(pfm)

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    nsites = len(sites)
    write_meme(output_dir, tag, pfm, background, nsites)
    write_pwm(output_dir, tag, pwm)
    write_pfm(output_dir, tag, pfm)
    write_sites(output=output_dir, tag=tag, sites=sites)
    with open(output_dir + '/fprs.txt', 'w') as file:
        for i in fprs:
            file.write(str(i) + '\n')
    file.close()
    return(0)
