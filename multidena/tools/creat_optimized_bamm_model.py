import subprocess
import os
import shutil
import random
from operator import itemgetter
from multidena.lib.common import read_peaks, write_fasta, read_bamm, \
creat_background, calculate_particial_auc, \
score_bamm, complement, make_pcm, make_pfm, write_meme, write_auc, \
write_auc_with_order, calculate_merged_roc, write_table, \
calculate_short_roc, write_roc, calculate_fprs, calculate_prc
from multidena.lib.speedup import creat_table_bootstrap


def run_streme(fasta_path, backgroud_path, dir_out, motif_length):
    args = ['streme', '--p', fasta_path,
            '--n', backgroud_path,
           '--objfun', 'de',
           '--minw',  str(motif_length),
           '--maxw',  str(motif_length),
           #'--w', str(motif_length),
           '-nmotifs', '5',
           '--text']
    p = subprocess.run(args, shell=False, capture_output=True)
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
    with open(f'{dir_out}/streme.txt', 'wb')  as file:
        file.write(p.stdout)
    return(0)


def run_streme_hmm_background(fasta_path, dir_out, motif_length):
    args = ['streme', '--p', fasta_path,
            '--kmer', '4',
           '--objfun', 'de',
           '--minw',  str(motif_length),
           '--maxw',  str(motif_length),
           #'--w', str(motif_length),
           '-nmotifs', '5',
           '--text']
    p = subprocess.run(args, shell=False, capture_output=True)
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
    with open(f'{dir_out}/streme.txt', 'wb')  as file:
        file.write(p.stdout)
    return(0)


def parse_streme(path):
    pfm = {'A': [], 'C': [], 'G': [], 'T': []}
    with open(path) as file:
        for line in file:
            if line.startswith('Background'):
                line = file.readline().strip().split()
                background = {'A': float(line[1]),
                             'C': float(line[3]),
                             'G': float(line[5]),
                             'T': float(line[7])}
            elif line.startswith('letter-probability'):
                break
        line = line.strip().split()
        length = int(line[5])
        nsites = int(line[7])
        for i in range(length):
            line = file.readline().strip().split()
            for letter, value in zip(pfm.keys(), line):
                    pfm[letter].append(float(value))
    file.close()
    return pfm, background, length, nsites


def create_bamm_model(peaks_path, backgroud_path, directory, order, meme, extend, basename):
    args = ['BaMMmotif', directory, peaks_path, '--PWMFile', meme,
            '--EM', '--order', str(order), '--Order', str(order),
           '--extend', str(extend), '--basename', str(basename),
           '--negSeqFile', backgroud_path]
    r = subprocess.run(args, capture_output=True)
    bamm_path = directory + '/{}_motif_1.ihbcp'.format(basename)
    bg_path = directory + '/{}.hbcp'.format(basename)
    bamm, order = read_bamm(bamm_path, bg_path)
    return(bamm, order)



def false_scores_bamm(peaks, bamm, order, length_of_site):
    false_scores = []
    append = false_scores.append
    for peak in peaks:
        complement_peak = complement(peak)
        full_peak = peak + 'N' * length_of_site + complement_peak
        n = len(full_peak) - length_of_site + 1
        for i in range(n):
            site = full_peak[i:length_of_site + i]
            if not len(set(site) - {'A', 'C', 'G', 'T'}) == 0:
                continue
            score = score_bamm(site, bamm, order, length_of_site)
            false_scores.append(score)
    return(false_scores)


def best_scores_bamm(peaks, bamm, order, length_of_site):
    true_scores = []
    for peak in peaks:
        complement_peak = complement(peak)
        best = -1000000
        full_peak = peak + 'N' * length_of_site + complement_peak
        n = len(full_peak) - length_of_site + 1
        for i in range(n):
            site = full_peak[i:length_of_site + i]
            if not len(set(site) - {'A', 'C', 'G', 'T'}) == 0:
                continue
            score = score_bamm(site, bamm, order, length_of_site)
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


def learn_optimized_bamm_support(peaks_path, backgroud_path, counter, order, length, pwm_auc_dir, tmp_dir, output_auc, pfpr):
    true_scores = []
    false_scores_roc = []
    false_scores_prc = []
    peaks = read_peaks(peaks_path)
    for step in ['odd', 'even']:
        meme = pwm_auc_dir + '/pwm_model_{0}_{1}.meme'.format(step, length)
        if step == 'odd':
            train_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 != 0]
            test_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 == 0]
        else:
            train_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 == 0]
            test_peaks = [p for index, p in enumerate(peaks, 1) if index % 2 != 0]
        write_fasta(train_peaks, tmp_dir + '/train.fasta')
        if os.path.isfile(backgroud_path):
            shuffled_peaks = read_peaks(backgroud_path)
            bamm, order = create_bamm_model(tmp_dir + '/train.fasta', backgroud_path, tmp_dir, order, meme, 0, length)
        else:
            shuffled_peaks = creat_background(test_peaks, length, counter)
            write_fasta(shuffled_peaks, tmp_dir + '/background.fasta')
            bamm, order = create_bamm_model(tmp_dir + '/train.fasta', tmp_dir + '/background.fasta', tmp_dir, order, meme, 0, length)
        for true_score in best_scores_bamm(test_peaks, bamm, order, length):
            true_scores.append(true_score)
        for false_score in false_scores_bamm(shuffled_peaks, bamm, order, length):
            false_scores_roc.append(false_score)
        for false_score in best_scores_bamm(shuffled_peaks, bamm, order, length):
            false_scores_prc.append(false_score)
        shutil.copy(tmp_dir + '/{}_motif_1.ihbcp'.format(length),
               output_auc + '/bamm_model_{0}_{1}_{2}.ihbcp'.format(step, order, length))
        shutil.copy(tmp_dir + '/{}.hbcp'.format(length),
               output_auc + '/bamm_{0}_{1}_{2}.hbcp'.format(step, order, length))
    fprs = calculate_fprs(true_scores, false_scores_roc)
    prc = calculate_prc(true_scores, false_scores_prc)
    roc = calculate_short_roc(fprs, step=1)
    merged_roc = calculate_merged_roc(fprs)
    auc_roc = calculate_particial_auc(merged_roc['TPR'], merged_roc['FPR'], pfpr)
    auc_prc = calculate_particial_auc(prc['PRECISION'], prc['RECALL'], 1.01)
    print("Length {0}; Order {1}".format(length, order), "pAUC at {0} = {1};".format(pfpr, auc_roc), "PRC = {0}".format(auc_prc))
    write_table(output_auc + '/statistics.txt', auc_roc, auc_prc, order, length)
    write_roc(output_auc + "/training_roc_merged_{0}_{1}.txt".format(length, order), merged_roc)
    write_roc(output_auc + "/training_roc_{0}_{1}.txt".format(length, order), roc)
    write_roc(output_auc + "/training_prc_{0}_{1}.txt".format(length, order), prc)
    return(0)


def learn_optimized_bamm(peaks_path, backgroud_path, counter, pwm_auc_dir, tmp_dir, output_auc, pfpr, min_length, max_length, length_step):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    if not os.path.isdir(output_auc):
        os.mkdir(output_auc)
    if os.path.exists(output_auc + '/statistics.txt'):
        os.remove(output_auc + '/statistics.txt')
    for order in range(1,5):
        #for length in range(12, 41, 4):
        for length in range(min_length, max_length + 1, length_step):
            learn_optimized_bamm_support(peaks_path, backgroud_path, counter, order, length, pwm_auc_dir, tmp_dir, output_auc, pfpr)
    pass


def choose_best_model(output_auc):
    auc = []
    with open(output_auc + '/statistics.txt') as file:
        for line in file:
            auc.append(tuple(map(float, line.strip().split())))
        file.close()
    auc.sort(key=itemgetter(-1))
    order = int(auc[-1][1])
    length = int(auc[-1][0])
    return(length, order)


def de_novo_with_oprimization_bamm(peaks_path, backgroud_path, pwm_auc_dir, tmp_dir,
    output_dir, output_auc, pfpr, min_length, max_length, length_step):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    else:
        shutil.rmtree(tmp_dir)
        os.mkdir(tmp_dir)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    if not os.path.isdir(output_auc):
        os.mkdir(output_auc)
    counter = 1000000
    learn_optimized_bamm(peaks_path, backgroud_path, counter, pwm_auc_dir, tmp_dir, output_auc, pfpr, min_length, max_length, length_step)
    length, order = choose_best_model(output_auc)
    meme = pwm_auc_dir + '/pwm_model_even_{}.meme'.format(length)
    if os.path.isfile(backgroud_path):
        run_streme(peaks_path, backgroud_path, tmp_dir, length)
        pfm, background_pfm, length_pfm, nsites = parse_streme(tmp_dir + '/streme.txt')
        tag = 'pfm_model'
        write_meme(tmp_dir, tag, pfm, background_pfm, nsites)
        meme = tmp_dir + '/pfm_model.meme'
        create_bamm_model(peaks_path, backgroud_path, tmp_dir, order, meme, 0, length)
    else:
        peaks = read_peaks(peaks_path)
        shuffled_peaks = creat_background(peaks, length, counter)
        write_fasta(shuffled_peaks, tmp_dir + '/background.fasta')
        run_streme_hmm_background(peaks_path, tmp_dir, length)
        pfm, background_pfm, length_pfm, nsites = parse_streme(tmp_dir + '/streme.txt')
        tag = 'pfm_model'
        write_meme(tmp_dir, tag, pfm, background_pfm, nsites)
        meme = tmp_dir + '/pfm_model.meme'
        create_bamm_model(peaks_path, tmp_dir + '/background.fasta', tmp_dir, order, meme, 0, length)
    shutil.copy(tmp_dir + '/pfm_model.meme',
           output_dir + '/pfm_model.meme')
    shutil.copy(tmp_dir + '/{}_motif_1.ihbcp'.format(length),
           output_dir + '/bamm_model.ihbcp')
    shutil.copy(tmp_dir + '/{}.hbcp'.format(length),
           output_dir + '/bamm.hbcp')
    shutil.copy(output_auc + '/training_prc_{0}_{1}.txt'.format(length, order),
             output_dir + '/prc.txt')
    shutil.copy(output_auc + '/training_roc_{0}_{1}.txt'.format(length, order),
             output_dir + '/roc.txt')
    shutil.copy(output_auc + '/training_roc_merged_{0}_{1}.txt'.format(length, order),
             output_dir + '/roc_merged.txt')
    shutil.rmtree(tmp_dir)
    return(length, order)
