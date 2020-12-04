import csv
from lib.common import read_fasta, read_dipwm
from lib.speedup import score_dipwm


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


def scan_seqs_by_dipwm(record, dipwm, threshold):
    results = []
    reverse_record = complement(record)
    length_dipwm = len(dipwm['AA']) + 1
    seq = record['seq']
    reverse_seq = reverse_record['seq']

    # first strand
    for i in range(len(seq) - length_dipwm + 1):
        site_seq = seq[i:length_dipwm + i]
        if not check_nucleotides(site_seq):
            continue
        s = score_dipwm(site_seq, dipwm)
        if s >= threshold:
            site_dict = dict()
            site_dict['name'] = record['name']
            site_dict['chromosome'] = record['chromosome']
            site_dict['start'] = str(int(record['start']) + i)
            site_dict['end'] = str(int(record['start']) + i + length_dipwm)
            site_dict['site'] = site_seq
            site_dict['strand'] = record['strand']
            site_dict['score'] = s
            results.append(site_dict)

    # second strand
    for i in range(len(seq) - length_dipwm + 1):
        site_seq = reverse_seq[i:length_dipwm + i]
        if not check_nucleotides(site_seq):
            continue
        s = score_dipwm(site_seq, dipwm)
        if s >= threshold:
            site_dict = dict()
            site_dict['name'] = record['name']
            site_dict['chromosome'] = record['chromosome']
            site_dict['start'] = str(int(record['end']) - i - length_dipwm)
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


def scan_by_dipwm(fasta_path, dipwm_path, threshold, results_path):
    fasta = read_fasta(fasta_path)
    dipwm = read_dipwm(dipwm_path)
    results = []
    for record in fasta:
      results += scan_seqs_by_dipwm(record, dipwm, threshold)
    write_csv(results_path, results)
    return(0)
