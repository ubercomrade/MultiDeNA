import re
import math
from operator import itemgetter

def read_bed(path_bed, fpr_thr_table, model):
    container = {}
    with open(path_bed) as file:
        for line in file:
            chrom, start, end, peak, score, strand, site = line.strip().split()
            if not peak in container:
                container[peak] = [{'chr': chrom,
                                    'start': int(start),
                                    'end': int(end),
                                    '-log10fpr': get_fpr(fpr_thr_table, float(score)),
                                    'strand': strand,
                                    'site': site.lower(),
                                    'model': model}]
            else:
                container[peak].append({'chr': chrom,
                                        'start': int(start),
                                        'end': int(end),
                                        '-log10fpr': get_fpr(fpr_thr_table, float(score)),
                                        'strand': strand,
                                        'site': site.lower(),
                                        'model': model})
    file.close()
    return(container)


def read_fasta(path):
    fasta = list()
    with open(path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                line = line[1:].strip().split(':')
                record = dict()
                record['name'] = line[0]
                record['chr'] = line[2]
                coordinates_strand = line[3]
                start, end = re.findall(r'\d*-\d*', coordinates_strand)[0].split('-')
                record['start'] = int(start)
                record['end'] = int(end)
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


def read_fpr_thr_table(path):
    container = list()
    append = container.append
    with open(path, 'r') as file:
        file.readline()
        for line in file:
            append(tuple(map(float, line.strip().split())))
    file.close()
    container = sorted(container, key=itemgetter(1))
    return(container)


def get_fpr(table, score):
    last_score, last_fpr = table[0]
    for line in table:
        if line[0] > score:
            break
        else:
            last_score, last_fpr = line
    return(-math.log10(last_fpr))


def write_mcot_format(path, bed, fasta, threshold):
    with open(path, 'w') as file:
        for index, line in enumerate(fasta):
            file.write('>peaks_{0}::{1}:{2}-{3}({4})\tSEQ {5}\tTHR {6}\n'.format(index, 
                                                                                 line['chr'],
                                                                                 line['start'],
                                                                                 line['end'],
                                                                                 line['strand'],
                                                                                 index + 1,
                                                                                 threshold))
            if line['name'] in bed:
                bed[line['name']] = sorted(bed[line['name']], key=itemgetter('start'))
                for site in bed[line['name']]:
                    pos = site['start'] - line['start']
                    file.write('{0}\t{1}\t{2}\t{3}\n'.format(pos, site['score'], site['strand'], site['site']))
    return(0)


def write_mcot_format(path, bed, fasta):
    with open(path, 'w') as file:
        for index, line in enumerate(fasta):
            file.write('>peaks_{0}::{1}:{2}-{3}({4})\tSEQ {5}\n'.format(index, 
                                                                                 line['chr'],
                                                                                 line['start'],
                                                                                 line['end'],
                                                                                 line['strand'],
                                                                                 index + 1))
            if line['name'] in bed:
                bed[line['name']] = sorted(bed[line['name']], key=itemgetter('start'))
                for site in bed[line['name']]:
                    start = site['start'] - line['start']
                    end = site['end'] - line['start']
                    file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(start, end, site['-log10fpr'],
                                                             site['strand'], site['site'],
                                                            site['model'], site['cluster']))
    return(0)


def combine_results(fasta_path, list_bed_path, list_path_fpr_table, list_models, path_to_write):
    bed = {}
    fasta = read_fasta(fasta_path)
    for model, bed_path, table_path in zip(list_models, list_bed_path, list_path_fpr_table):
        table = read_fpr_thr_table(table_path)
        for key, value in read_bed(bed_path, table, model).items():
            if key in bed:
                bed[key].extend(value)
            else:
                bed[key] = value
    for k in bed.keys():
        bed[k] = sorted(bed[k], key=itemgetter('start'))

    for k in bed.keys():
        peak = bed[k]
        index = 1
        site = peak[0]
        site['cluster'] = index
        models = [site['model']]
        for i in peak[1:]:
            if (i['start'] < site['end']) and (i['end'] > site['start']) and (not i['model'] in models):
                i['cluster'] = index
                models.append(i['model'])
            else:
                index += 1
                i['cluster'] = index
                models = [i['model']]
            site = i
    write_mcot_format(path_to_write, bed, fasta)
    return(0)


