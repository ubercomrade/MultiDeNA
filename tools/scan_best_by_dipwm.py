import re
from lib.speedup import score_dipwm
from lib.common import read_dipwm


def read_fasta(path):
    fasta = list()
    with open(path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                line = line[1:].strip().split(':')
                record = dict()
                record['name'] = line[0]
                record['chromosome'] = line[2]
                coordinates_strand = line[3]

                start, end = re.findall(r'\d*-\d*', coordinates_strand)[0].split('-')
                record['start'] = start
                record['end'] = end

                strand = re.findall(r'\(.\)', coordinates_strand[:-3])
                if not strand == []:
                    record['strand'] = strand[0].strip('()')
                else:
                    record['strand'] = '+'
            else:
                record['seq'] = line.strip().upper()
                fasta.append(record)
    file.close()
    return(fasta)


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


def scan_seq_by_dipwm(record, dipwm):
    results = []
    reverse_record = complement(record)
    length_dipwm = len(dipwm['AA']) + 1
    seq = record['seq']
    reverse_seq = reverse_record['seq']
    threshold = -1000000

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
            threshold = s

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
            threshold = s

    results.append(site_dict)
    return(results)


def write_list(path, data):
    scores = [i['score'] for i in data]
    with open(path, "w") as file:
        for line in scores:
            file.write("{0}\n".format(line))
    file.close()
    pass


def scan_best_by_dipwm(results_path, dipwm_path, fasta_path):
    fasta = read_fasta(fasta_path)
    dipwm = read_dipwm(dipwm_path)
    results = []
    for record in fasta:
      results += scan_seq_by_dipwm(record, dipwm)
    write_list(results_path, results)
    return(0)
