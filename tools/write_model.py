import math
import itertools
from collections import Counter


def read_sites(path):
    with open(path, 'r') as file:
        sequences = [i.strip().upper() for i in file if i.strip()[0] != '>']
        sequences = [i for i in sequences if not 'N' in i]
    return(sequences)


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


def make_pfm_from_pcm(pcm):
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


def make_pwm_from_pcm(pcm, background):
    pwm = {}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        pwm[i[0]] = []

    pfm = make_pfm_from_pcm(pcm)
    first_key = list(pcm.keys())[0]
    for i in range(len(pfm[first_key])):
        for j in pfm.keys():
            pwm[j].append(math.log2(pfm[j][i] / background[j]))
    return(pwm)


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


def check_legth_of_sites(seq):
    length = [len(i) for i in seq]
    c = Counter(length).most_common(1)[0][0]
    seq = [i for i in seq if len(i) == c]
    return(seq)


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


def write_model(fasta_path, output_dir, tag):
    seq = read_sites(fasta_path)
    if len(seq) > 1:
        seq = list(set(seq))
        seq = check_legth_of_sites(seq)
        background = {'A': 0.25,
                     'C': 0.25,
                     'G': 0.25,
                     'T': 0.25}
        nsites = len(seq)

        pcm = make_pcm(seq)
        pfm = make_pfm_from_pcm(pcm)
        pwm = make_pwm_from_pcm(pcm, background)

        write_pwm(output_dir, tag, pwm)
        write_pfm(output_dir, tag, pfm)
        write_meme(output_dir, tag, pfm, background, nsites)
    else:
        print("INMODE model didn't fined sites -> cant't creat matrix models from sites")
    return(0)