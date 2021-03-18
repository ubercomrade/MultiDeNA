import csv
import pickle
import numpy as np
from strum import strum
from lib.common import read_fasta, read_strum


def complement(record):
    output = dict(record)
    strand = record['strand']
    seq = str()
    if strand == '+':
        output['strand'] = '-'
    else:
        output['strand'] = '+'
    seq = output['seq'].replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]
    output['seq'] = seq
    return(output)


def check_nucleotides(site):
    s = set(site)
    n = {'A', 'C', 'G', 'T'}
    if len(s - n) == 0:
        return(True)
    else:
        return(False)


def scan_seqs_by_strum(record, strum_model, length, threshold):
    results = []
    reverse_record = complement(record)
    seq = record['seq']
    reverse_seq = reverse_record['seq']

    # first strand
    scores = strum_model.score_seq(seq)
    upper_threshold = np.nonzero(scores >= threshold)[0]
    for i in upper_threshold:
        s = scores[i]
        site_seq = seq[i:length + i]
        site_dict = dict()
        site_dict['name'] = record['name']
        site_dict['chromosome'] = record['chromosome']
        site_dict['start'] = str(int(record['start']) + i)
        site_dict['end'] = str(int(record['start']) + i + length)
        site_dict['site'] = site_seq
        site_dict['strand'] = record['strand']
        site_dict['score'] = s
        results.append(site_dict)

    # second strand
    scores = strum_model.score_seq(reverse_seq)
    upper_threshold = np.nonzero(scores >= threshold)[0]
    for i in upper_threshold:
        s = scores[i]
        site_seq = reverse_seq[i:length + i]
        site_dict = dict()
        site_dict['name'] = record['name']
        site_dict['chromosome'] = record['chromosome']
        site_dict['start'] = str(int(record['end']) - i - length)
        site_dict['end'] = str(int(record['end']) - i)
        site_dict['site'] = site_seq
        site_dict['strand'] = reverse_record['strand']
        site_dict['score'] = s
        results.append(site_dict)
    return(results)


def write_csv(path, data):
    with open(path, 'w') as csvfile:
        fieldnames = ['chromosome', 'start', 'end', 'name', 'score', 'strand', 'site']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        for line in data:
            writer.writerow(line)
    pass


def scan_by_strum(fasta_path, strum_path, threshold, results_path):
    fasta = read_fasta(fasta_path)
    strum_model = read_strum(strum_path)
    length = strum_model.k
    results = []
    for record in fasta:
        results += scan_seqs_by_strum(record, strum_model, length, threshold)
    write_csv(results_path, results)
    return(0)
