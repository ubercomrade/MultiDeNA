import random
import shutil
import os
import subprocess
from lib.common import read_peaks, sites_to_dipwm, \
creat_background, write_fasta, complement, \
write_roc, calculate_fprs, calculate_short_roc, calculate_merged_roc
from lib.speedup import creat_table_bootstrap, score_dipwm


def run_chipmunk(path_to_java, path_to_chipmunk, fasta_path, path_out, motif_length_start, motif_length_end, cpu_count):
    args = [path_to_java, '-cp', path_to_chipmunk,
                   'ru.autosome.di.ChIPMunk', str(motif_length_start), str(motif_length_end), 'yes', '1.0',
                   's:{}'.format(fasta_path),
                  '100', '10', '1', str(cpu_count), 'random']
    p = subprocess.run(args, shell=False, capture_output=True)
    out = p.stdout
    with open(path_out, 'wb') as file:
        file.write(out)
    return(0)


def parse_chipmunk(path):
    with open(path, 'r') as file:
        container = []
        for line in file:
            d = {'name': str(), 'start': int(), 'end': int(),
                 'seq': str(), 'strand': str()}
            if line.startswith('WORD|'):
                line = line[5:].strip()
                line = line.split()
                d['name'] = 'peaks_' + str(int(line[0]) - 1)
                d['start'] = int(line[1])
                d['end'] = int(line[1]) + len(line[2])
                d['seq'] = line[2]
                if line[4] == 'direct':
                    d['strand'] = '+'
                else:
                    d['strand'] = '-'
                container.append(d)
            else:
                continue
    seqs = [i['seq'] for i in container if not 'N' in i['seq']]
    return(seqs)


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


def bootstrap_dipwm(peaks, length_of_site, counter, path_to_java, path_to_chipmunk, tmp_dir, cpu_count):
    number_of_peaks = len(peaks)
    fpr_of_every_site = []
    for i in range(5):
        true_scores = []
        false_scores = []
        train_peaks = random.sample(peaks, k=round(0.9 * number_of_peaks))
        test_peaks = [peak for peak in peaks if not peak in train_peaks]
        shuffled_peaks = creat_background(test_peaks, length_of_site, counter)
        write_fasta(train_peaks, tmp_dir + '/train.fasta')
        run_chipmunk(path_to_java, path_to_chipmunk,
                     tmp_dir + '/train.fasta', tmp_dir + '/chipmunk_results.txt',
                     length_of_site, length_of_site, cpu_count)
        sites = parse_chipmunk(tmp_dir + '/chipmunk_results.txt')
        sites = list(set(sites))
        dipwm = sites_to_dipwm(sites)
        for true_score in true_scores_dipwm(test_peaks, dipwm, length_of_site):
            true_scores.append(true_score)
        for false_score in false_scores_dipwm(shuffled_peaks, dipwm, length_of_site):
            false_scores.append(false_score)
        fpr_of_every_site += calculate_fprs(true_scores, false_scores)
    fpr_of_every_site.sort()
    return(fpr_of_every_site)



def bootstrap_for_dipwm(peaks_path, results_path, results_path_wide, length_of_site, path_to_java, path_to_chipmunk, tmp_dir, cpu_count, counter=5000000):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    peaks = read_peaks(peaks_path)
    fprs = bootstrap_dipwm(peaks, length_of_site, counter, path_to_java, path_to_chipmunk, tmp_dir, cpu_count)
    short_roc = calculate_short_roc(fprs, step=1)
    merged_roc = calculate_merged_roc(fprs)
    shutil.rmtree(tmp_dir)
    write_roc(results_path, short_roc)
    write_roc(results_path_wide, merged_roc)
    return(0)