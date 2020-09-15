import subprocess
import os
import shutil
import random
from lib.common import read_peaks, write_fasta, read_bamm, \
write_table_bootstrap, creat_background, \
score_bamm, complement, make_pcm, make_pfm, write_meme, \
write_table_bootstrap_wide
from lib.speedup import creat_table_bootstrap


def create_bamm_model(directory, order, meme):
    fasta_path = directory + '/train.fasta'
    args = ['BaMMmotif', directory, fasta_path, '--PWMFile', meme, '--EM', '--order', str(order), '--Order', str(order)]
    r = subprocess.run(args, capture_output=True)
    bamm_path = directory + '/train_motif_1.ihbcp'
    bg_path = directory + '/train.hbcp'
    bamm, order = read_bamm(bamm_path, bg_path)
    return(bamm, order)


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



def false_scores_bamm(peaks, bamm, order, length_of_site):
    false_scores = []
    append = false_scores.append
    for peak in peaks:
        complement_peak = complement(peak)
        full_peak = peak + 'N' * length_of_site + complement_peak
        n = len(full_peak) - length_of_site + 1
        for i in range(n):
            site = full_peak[i:length_of_site + i]
            if 'N' in site:
                continue
            score = score_bamm(site, bamm, order, length_of_site)
            false_scores.append(score)
    return(false_scores)


def true_scores_bamm(peaks, bamm, order, length_of_site):
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
            score = score_bamm(site, bamm, order, length_of_site)
            if score >= best:
                best = score
        true_scores.append(best)
    return(true_scores)


def bootstrap_bamm(peaks, length_of_site, counter, order, path_to_chipmunk, path_to_java, cpu_count, tmp_dir):
    true_scores = []
    false_scores = []
    number_of_peaks = len(peaks)
    background = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}
    for i in range(5):
        train_peaks = random.choices(peaks, k=round(0.9 * number_of_peaks))
        test_peaks = [peak for peak in peaks if not peak in train_peaks]
        shuffled_peaks = creat_background(test_peaks, length_of_site, counter / 5)
        write_fasta(train_peaks, tmp_dir + '/train.fasta')
        run_chipmunk(path_to_java, path_to_chipmunk,
                    tmp_dir + '/train.fasta', tmp_dir + '/chipmunk_results.txt',
                    length_of_site, length_of_site, cpu_count)
        sites = parse_chipmunk(tmp_dir + '/chipmunk_results.txt')
        sites = list(set(sites))
        nsites = len(sites)
        pcm = make_pcm(sites)
        pfm = make_pfm(pcm)
        write_meme(tmp_dir, "anchor", pfm, background, nsites)
        bamm, order = create_bamm_model(tmp_dir, order, tmp_dir + '/anchor.meme')
        for true_score in true_scores_bamm(test_peaks, bamm, order, length_of_site):
            true_scores.append(true_score)
        for false_score in false_scores_bamm(shuffled_peaks, bamm, order, length_of_site):
            false_scores.append(false_score)
    table = creat_table_bootstrap(true_scores, false_scores)
    table_full = calculate_roc(true_scores, false_scores)
    return(table, table_full)


def bootstrap_for_bamm(peaks_path, results_path, results_path_wide, length_of_site, 
                       path_to_chipmunk, path_to_java, cpu_count, 
                       tmp_dir, counter = 5000000, order=2):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    peaks = read_peaks(peaks_path)
    table, table_full = bootstrap_bamm(peaks, length_of_site, counter, order, path_to_chipmunk, path_to_java, cpu_count, tmp_dir)
    write_table_bootstrap(results_path, table)
    write_table_bootstrap_wide(results_path_wide, table_full)
    shutil.rmtree(tmp_dir)
    return(0)



