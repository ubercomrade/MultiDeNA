import subprocess
import os
import random
import math
import shutil
from operator import itemgetter
from lib.common import read_peaks, write_table_bootstrap, \
creat_background, complement, write_table_bootstrap_wide, \
calculate_roc
from lib.speedup import creat_table_bootstrap


def read_inmode_bed(path):
    table = []
    with open(path) as file:
        for line in file:
            line = line.strip().split('\t')
            line[0] = int(line[0])
            line[1] = int(line[1])
            line[2] = int(line[2])
            line[4] = float(line[4])
            table.append(line)
    return(table)


def write_fasta(sites, tmp_dir, tag):
    with open('{0}/{1}.fa'.format(tmp_dir, tag), 'w') as file:
        for index, site in enumerate(sites):
            file.write('>{0}\n{1}\n'.format(index, site))
    return(0)


def true_scores_inmode(path_to_inmode, path_to_java, motif_length, order, tmp_dir, tag):
    scores = []
    args = [path_to_java, '-Xmx16G', '-Xms1G', 
            '-jar',
            path_to_inmode, 'scan',
            'i={0}/Learned_DeNovo({1},{2},2)_motif/XML_of_DeNovo({1},{2},2)_motif.xml'.format(tmp_dir, motif_length, order),
            'id={0}/{1}.fa'.format(tmp_dir, tag), 'f=1.0', 'outdir={}'.format(tmp_dir), 'bs=false']
    r = subprocess.run(args, capture_output=True)
    scores = []
    table = read_inmode_bed('{0}/{1}'.format(tmp_dir, "/Motif_hits_from_SequenceScan(1.0).BED"))
    table.sort(key=itemgetter(0, 4))
    last_index = 0
    for line in table:
        index = line[0]
        score = line[4]
        if last_index != index:
            scores.append(last_score)
        last_score = score
        last_index = index
    scores.append(score)
    scores = [math.log(float(i), 10) for i in scores]
    os.remove(tmp_dir + '/Binding_sites_from_SequenceScan(1.0).txt')
    os.remove(tmp_dir + '/Motif_hits_from_SequenceScan(1.0).BED')
    os.remove(tmp_dir + '/protocol_scan.txt')
    os.remove(tmp_dir + '/{}.fa'.format(tag))
    return(scores)


def false_scores_inmode(path_to_inmode, path_to_java, motif_length, order, tmp_dir, tag):
    scores = []
    args = [path_to_java, '-Xmx16G', '-Xms1G', 
            '-jar',
            path_to_inmode, 'scan',
            'i={0}/Learned_DeNovo({1},{2},2)_motif/XML_of_DeNovo({1},{2},2)_motif.xml'.format(tmp_dir, motif_length, order),
            'id={0}/{1}.fa'.format(tmp_dir, tag), 'f=1.0', 'outdir={}'.format(tmp_dir), 'bs=false']
    r = subprocess.run(args, capture_output=True)
    with open('{0}/{1}'.format(tmp_dir, "/Motif_hits_from_SequenceScan(1.0).BED")) as file:
        for line in file:
            scores.append(math.log(float(line.split()[4]), 10))
    os.remove(tmp_dir + '/Binding_sites_from_SequenceScan(1.0).txt')
    os.remove(tmp_dir + '/Motif_hits_from_SequenceScan(1.0).BED')
    os.remove(tmp_dir + '/protocol_scan.txt')
    os.remove(tmp_dir + '/{}.fa'.format(tag))
    return(scores)


def make_inmode(path_to_inmode, path_to_java, motif_length, order, tmp_dir):
    args = [path_to_java, '-Xmx16G', '-Xms1G', '-jar', path_to_inmode,
    'denovo', 'i={}/train.fa'.format(tmp_dir), 'm={}'.format(motif_length), 'outdir={}'.format(tmp_dir),
    'mo={}'.format(order)]
    r = subprocess.run(args, capture_output=True)
    return(0)


def bootstrap_inmode(peaks, length_of_site, counter, path_to_inmode, path_to_java, tmp_dir, order):
    true_scores = []
    false_scores = []
    number_of_peaks = len(peaks)
    for i in range(5):
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)
        train_peaks = random.choices(peaks, k=round(0.9 * number_of_peaks))
        test_peaks = [peak for peak in peaks  if not peak in train_peaks]
        shuffled_peaks = creat_background(test_peaks, length_of_site, counter / 5)
        write_fasta(train_peaks, tmp_dir, "train")
        write_fasta(test_peaks, tmp_dir, "test")
        write_fasta(shuffled_peaks, tmp_dir, "shuffled")
        make_inmode(path_to_inmode, path_to_java, length_of_site, order, tmp_dir)
        for true_score in true_scores_inmode(path_to_inmode, path_to_java, length_of_site, order, tmp_dir, "test"):
            true_scores.append(true_score)
        for false_score in false_scores_inmode(path_to_inmode, path_to_java, length_of_site, order, tmp_dir, "shuffled"):
            false_scores.append(false_score)
        shutil.rmtree(tmp_dir)
    table = creat_table_bootstrap(true_scores, false_scores)
    table_full = calculate_roc(true_scores, false_scores)
    return(table, table_full)


def bootstrap_for_inmode(peaks_path, results_path, results_path_wide, length_of_site, path_to_inmode, path_to_java, tmp_dir, counter=5000000, order=2):
    peaks = read_peaks(peaks_path)
    table, table_full = bootstrap_inmode(peaks, length_of_site, counter, path_to_inmode, path_to_java, tmp_dir, order)
    write_table_bootstrap(results_path, table)
    write_table_bootstrap_wide(results_path_wide, table_full)
    return(0)

