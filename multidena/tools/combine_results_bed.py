import re
import math
import copy
import random
from operator import itemgetter
from multidena.lib.speedup import score_pwm
from multidena.lib.common import make_pcm, make_pfm, make_pwm



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
        if line[0] < score:
            break
        else:
            last_score, last_fpr = line
    return(-math.log10(last_fpr))


## bed_combine_format
def sites_to_pwm(peaks, model):
    sites = []
    for k in peaks.keys():
        peak = peaks[k]
        for i in peak:
            if i['model'] == model:
                sites.append(i['site'].upper())
    pcm = make_pcm(sites)
    pfm = make_pfm(pcm)
    pwm = make_pwm(pfm)
    return(pwm, sites)


def score_pwm(seq, pwm):
    score = 0
    for position, letter in enumerate(seq):
        score += pwm[letter][position]
    return(score)

def complement(seq):
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


def get_best_site_position_by_pwm(seq, pwm):
    results = []
    complement_seq = complement(seq)
    length = len(pwm['A'])
    threshold = -1000000
    # first strand
    for i in range(len(seq) - length + 1):
        site = seq[i:length + i]
        s = score_pwm(site, pwm)
        if s >= threshold:
            threshold = s
            position = i
            strand = '+'
            main_site = site
    # second strand
    for i in range(len(complement_seq) - length + 1):
        site = complement_seq[i:length + i]
        s = score_pwm(site, pwm)
        if s >= threshold:
            threshold = s
            position = i
            strand = '-'
            main_site = site
    return(position, main_site, length, strand)


def get_aligned_sites_by_pwm(peaks, pwm):
    peaks_copy = copy.deepcopy(peaks)
    for index, k in enumerate(peaks_copy.keys()):
        peak = peaks_copy[k]
        for i in peak:
            seq = i['site'].upper()
            position, main_site, length, strand = get_best_site_position_by_pwm(seq, pwm)
            i['site'] = main_site
            i['start'] += position
            i['end'] = i['start'] + length
            if strand == '-':
                if i['strand'] == '+':
                    i['strand'] = '-'
                else:
                    i['strand'] = '+'
    return(peaks_copy)


def calculate_aligned_positions(peaks, pwm):
    container = []
    for k in peaks.keys():
        peak = peaks[k]
        for i in peak:
            seq = i['site'].upper()
            position, main_site, length, strand = get_best_site_position_by_pwm(seq, pwm)
            i['align'] = position
            container.append(position)
            if strand == '-':
                i['site'] = complement(i['site'])
                if i['strand'] == '+':
                    i['strand'] = '-'
                else:
                    i['strand'] = '+'
    return(peaks, container)


def combine_results(fasta_path, list_bed_path, list_path_fpr_table, list_models):
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
    return(bed)


def merge_sites_in_clusters(bed):
    container = dict()
    for k in bed.keys():
        peak = bed[k]
        container[k] = []
        cluster = 1
        site = peak[0]
        record = {'start': site['start'],
                  'end': site['end'],
                  'site': site['site'],
                  'chr': site['chr'],
                  'models': [site['model']],
                  'scores': [site['-log10fpr']],
                  'strand': site['strand']}
        for site in peak[1:]:
            if site['cluster'] == cluster:
                record['models'].append(site['model'])
                record['scores'].append(site['-log10fpr'])
                if record['end'] < site['end']:
                    record['site'] = record['site'] + site['site'][record['end']-site['end']:]
                record['end'] = site['end']
            else:
                container[k].append(record)
                cluster += 1
                record = {'start': site['start'],
                          'end': site['end'],
                          'site': site['site'],
                          'chr': site['chr'],
                          'models': [site['model']],
                          'scores': [site['-log10fpr']],
                          'strand': site['strand']}
        container[k].append(record)
    return(container)


def choose_max_length_of_site(data):
    lengths = []
    for k in data.keys():
        peak = data[k]
        for i in peak:
            lengths.append(len(i['site']))
    lengths_counter = Counter(lengths)
    choosen_length = max([l for l in lengths_counter.keys() if lengths_counter[l] > 5])
    return(choosen_length)


def extend_sites(bed, fasta, max_length):
    bed_copy = bed.copy()
    for k in bed_copy.keys():
        for peak in bed_copy[k]:
            length = len(peak['site'])
            if length <= max_length:
                number_of_nulceotides_left = (max_length - length) // 2
                number_of_nulceotides_right = (max_length - length) - number_of_nulceotides_left
                number_of_peak = int(k.split('_')[1])
                relative_start = peak['start'] - fasta[number_of_peak]['start']
                relative_end = peak['end'] - fasta[number_of_peak]['start']
                if relative_start - number_of_nulceotides_left < 0:
                    number_of_nulceotides_left_corrected = number_of_nulceotides_left + (relative_start - number_of_nulceotides_left)
                    left_tail = ''.join(random.choices('ACGT', k=number_of_nulceotides_left - number_of_nulceotides_left_corrected)) \
                    + fasta[number_of_peak]['seq'][relative_start - number_of_nulceotides_left_corrected:relative_start].upper()
                    right_tail = fasta[number_of_peak]['seq'][relative_end:relative_end + number_of_nulceotides_right].upper()
                elif relative_end + number_of_nulceotides_right > fasta[number_of_peak]['end'] - fasta[number_of_peak]['start']:
                    left_tail = fasta[number_of_peak]['seq'][relative_start - number_of_nulceotides_left:relative_start].upper()
                    right_tail = fasta[number_of_peak]['seq'][relative_end:].upper() + \
                    ''.join(random.choices('ACGT', k=relative_end + number_of_nulceotides_right - (fasta[number_of_peak]['end'] - fasta[number_of_peak]['start'])))
                else:
                    left_tail = fasta[number_of_peak]['seq'][relative_start - number_of_nulceotides_left:relative_start].upper()
                    right_tail = fasta[number_of_peak]['seq'][relative_end:relative_end + number_of_nulceotides_right].upper()
                peak['site'] = left_tail + peak['site'] + right_tail
                peak['start'] = peak['start'] - number_of_nulceotides_left
                peak['end'] = peak['end'] + number_of_nulceotides_right
    return(bed_copy)


def get_flangs_for_sites(bed, fasta):
    bed_copy = bed.copy()
    number_of_nulceotides = 30
    for k in bed_copy.keys():
        number_of_peak = int(k.split('_')[1])
        for peak in bed_copy[k]:
            strand = peak['strand']
            relative_start = peak['start'] - fasta[number_of_peak]['start']
            relative_end = peak['end'] - fasta[number_of_peak]['start']
            left_tail = fasta[number_of_peak]['seq'][relative_start - number_of_nulceotides:relative_start].upper()
            right_tail = fasta[number_of_peak]['seq'][relative_end:relative_end + number_of_nulceotides].upper()
            if strand == '-':
                right_tail = complement(right_tail)
                left_tail = complement(left_tail)
            peak['right_tail'] = right_tail
            peak['left_tail'] = left_tail
    return(bed_copy)


def get_best_model(bed):
    bed_copy = bed.copy()
    for k in bed_copy.keys():
        for peak in bed_copy[k]:
            models = peak['models']
            scores = peak['scores']
            index = scores.index(max(scores))
            peak['best_model'] = models[index]
    return(bed_copy)


def bed_to_list(bed, max_length):
    container = []
    bed_copy = bed.copy()
    for k in bed_copy.keys():
        for peak in bed_copy[k]:
            peak['peak'] = k
            peak['models'] = ';'.join(peak['models'])
            peak['scores'] = ';'.join(map(str, peak['scores']))
            if len(peak['site']) <= max_length:
                container.append(peak)
    return(container)


def calculate_gc_at_content(bed):
    for line in bed:
        flangs = line['right_tail'] + line['left_tail']
        length = len(flangs)
        content = 0
        for n in flangs:
            if n == 'C' or n == 'G':
                content += 1
        line['GC'] = round(content / length, 2)
        line['AT'] = round(1 - line['GC'], 2)
    return(bed)


def write_wide_bed(bed, path):
    with open(path, 'w', newline='') as csvfile:
        #fieldnames = ['chr', 'start', 'end', 'peak', 'scores', 'strand', 'site', 'models', 'best_model', 'left_tail', 'right_tail', 'GC', 'AT']
        fieldnames = ['chr', 'start', 'end', 'peak', 'scores', 'strand', 'site', 'models', 'best_model']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for line in bed:
            writer.writerow(line)
    pass


def combine_results_bed_format(fasta_path, list_bed_path, list_path_fpr_table, 
    list_models, path_to_write_clusters, path_to_write_aligned):
    fasta = read_fasta(fasta_path)
    bed = combine_results(fasta_path, list_bed_path, list_path_fpr_table, list_models)
    pwm, sites = sites_to_pwm(bed, 'pwm')
    max_length = choose_max_length_of_site(bed)
    bed_clusters = merge_sites_in_clusters(bed)
    # bed, al_positions = calculate_aligned_positions(bed, pwm)
    # max_length = choose_max_length_of_site(bed)
    bed = get_aligned_sites_by_pwm(bed_clusters, pwm)
    bed = extend_sites(bed, fasta, max_length)
    # bed = get_flangs_for_sites(bed, fasta)
    # bed = calculate_gc_at_content(bed)
    bed = get_best_model(bed)
    #bed_clusters = bed_to_list(bed_clusters, max_length)
    bed = bed_to_list(bed, max_length)
    write_wide_bed(bed_clusters, path_to_write_clusters)
    write_wide_bed(bed, path_to_write_aligned)
    pass
