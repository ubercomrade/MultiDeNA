import random
import shutil
import os
import subprocess
from lib.common import read_peaks, sites_to_pwm, creat_background, write_fasta, complement, make_pcm, make_pfm, make_pwm, write_meme, write_pwm, write_pfm
from lib.speedup import creat_table_bootstrap, score_pwm


def run_chipmunk(path_to_java, path_to_chipmunk, fasta_path, path_out, motif_length_start, motif_length_end, cpu_count):
    args = [path_to_java, '-cp', path_to_chipmunk,
                   'ru.autosome.ChIPMunk', str(motif_length_start), str(motif_length_end), 'yes', '1.0',
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


def false_scores_pwm(peaks, pwm, length_of_site):
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
            score = score_pwm(site, pwm)
            false_scores.append(score)
    return(false_scores)


def true_scores_pwm(peaks, pwm, length_of_site):
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
            score = score_pwm(site, pwm)
            if score >= best:
                best = score
        true_scores.append(best)
    return(true_scores)


def fpr_at_tpr(true_scores, false_scores, tpr):
    true_scores.sort(reverse=True)
    false_scores.sort(reverse=True)
    false_length = len(false_scores)
    true_length = len(true_scores)
    score = true_scores[round(true_length * tpr) - 1]
    actual_tpr = sum([1 if true_score >= score else 0 for true_score in true_scores]) / true_length
    fpr = sum([1 if false_score >= score else 0 for false_score in false_scores]) / false_length
    return(fpr)


def write_sites(output, tag, sites):
    with open(output + '/' + tag + '.fasta', 'w') as file:
        for index, site in enumerate(sites):
            file.write(site + '\n')
    return(0)


def learn_optimized_pwm(peaks_path, counter, path_to_java, path_to_chipmunk, tmp_dir, cpu_count):

    # INITIAL STEP
    length = 12
    tpr = 0.3
    true_scores = []
    false_scores = []
    peaks = read_peaks(peaks_path)
    shuffled_peaks = creat_background(peaks, length, counter)
    run_chipmunk(path_to_java, path_to_chipmunk,
                 peaks_path, tmp_dir + '/chipmunk_results.txt',
                 length, length, cpu_count)
    sites_current = parse_chipmunk(tmp_dir + '/chipmunk_results.txt')
    sites_current = list(set(sites_current))
    pwm = sites_to_pwm(sites_current)
    for true_score in true_scores_pwm(peaks, pwm, length):
        true_scores.append(true_score)
    for false_score in false_scores_pwm(shuffled_peaks, pwm, length):
        false_scores.append(false_score)
    fpr_current = fpr_at_tpr(true_scores, false_scores, tpr)
    print(length, fpr_current)
    for length in range(14, 34, 2):
        true_scores = []
        false_scores = []
        peaks = read_peaks(peaks_path)
        shuffled_peaks = creat_background(peaks, length, counter)
        run_chipmunk(path_to_java, path_to_chipmunk,
                     peaks_path, tmp_dir + '/chipmunk_results.txt',
                     length, length, cpu_count)
        sites_new = parse_chipmunk(tmp_dir + '/chipmunk_results.txt')
        sites_new = list(set(sites_new))
        pwm = sites_to_pwm(sites_new)
        for true_score in true_scores_pwm(peaks, pwm, length):
            true_scores.append(true_score)
        for false_score in false_scores_pwm(shuffled_peaks, pwm, length):
            false_scores.append(false_score)
        fpr_new = fpr_at_tpr(true_scores, false_scores, tpr)
        if fpr_new < fpr_current: and (1 - fpr_new/fpr_current) * 100 > 5:
            sites_current = sites_new[:]
            fpr_current = fpr_new
            print(length, fpr_current)
        else:
            print(length, fpr_new)
            break
    return(sites_current, length, fpr_current)


def de_novo_with_oprimization_pwm(peaks_path, path_to_java, path_to_chipmunk, tmp_dir, output_dir, cpu_count):
    counter = 4000000
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    sites, length, fpr_current = learn_optimized_pwm(peaks_path, counter, path_to_java, path_to_chipmunk, tmp_dir, cpu_count)
    shutil.rmtree(tmp_dir)
    pcm = make_pcm(sites)
    pfm = make_pfm(pcm)
    pwm = make_pwm(pfm)

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    nsites = len(sites)
    background = {'A': 0.25,
                 'C': 0.25,
                 'G': 0.25,
                 'T': 0.25}
    tag = 'pwm_model'
    write_meme(output_dir, tag, pfm, background, nsites)
    write_pwm(output_dir, tag, pwm)
    write_pfm(output_dir, tag, pfm)
    write_sites(output=output_dir, tag=tag, sites=sites)
    return(0)
