import subprocess
import os
import random
import math
import shutil
from lib.common import read_seqs, \
write_table_bootstrap, \
creat_random_sample, \
creat_table_bootstrap


def write_fasta(sites, tmp_dir, tag):
    with open('{0}/{1}.fa'.format(tmp_dir, tag), 'w') as file:
        for index, site in enumerate(sites):
            file.write('>{0}\n{1}\n'.format(index, site))
    return(0)


def calculate_scores_inmode_bootsrap(path_to_inmode, path_to_java, motif_length, tmp_dir, tag):
    container = []
    args = [path_to_java, '-Xmx4096m', '-Xms1024m', 
            '-jar',
            path_to_inmode, 'scan',
            'i={0}/Learned_DeNovo({1},2,2)_motif/XML_of_DeNovo({1},2,2)_motif.xml'.format(tmp_dir, motif_length),
            'id={}/{}.fa'.format(tmp_dir, tag), 'f=1.0', 'outdir={}'.format(tmp_dir), 'bs=false']
    r = subprocess.run(args, capture_output=True)
    with open('{0}/{1}'.format(tmp_dir, "/Motif_hits_from_SequenceScan(1.0).BED")) as file:
        for line in file:
            container.append(math.log(float(line.split()[4]), 10))
    os.remove(tmp_dir + '/Binding_sites_from_SequenceScan(1.0).txt')
    os.remove(tmp_dir + '/Motif_hits_from_SequenceScan(1.0).BED')
    os.remove(tmp_dir + '/protocol_scan.txt')
    os.remove(tmp_dir + '/{}.fa'.format(tag))
    return(container)


def make_inmode(path_to_inmode, path_to_java, motif_length, order, tmp_dir):
    args = [path_to_java, '-Xmx4096m', '-Xms1024m', '-jar', path_to_inmode,
    'denovo', 'i={}/train.fa'.format(tmp_dir), 'm={}'.format(motif_length), 'outdir={}'.format(tmp_dir),
    'mo={}'.format(order)]
    r = subprocess.run(args, capture_output=True)
    return(0)


def bootstrap_inmode(sites, path_to_inmode, path_to_java, tmp_dir, order, size_of_random_sample):
    true_scores = []
    false_scores = []
    number_of_sites = len(sites)
    len_of_site = len(sites[0])

    for i in range(10):
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)
        train_sample = random.choices(sites, k=round(0.9 * number_of_sites))
        test_sample = [site for site in sites  if not site in train_sample]
        random_sample = creat_random_sample(test_sample, size_of_random_sample)
        write_fasta(train_sample, tmp_dir, "train")
        write_fasta(test_sample, tmp_dir, "test")
        write_fasta(random_sample, tmp_dir, "shuffled")
        make_inmode(path_to_inmode, path_to_java, len_of_site, order, tmp_dir)
        for true_score in calculate_scores_inmode_bootsrap(path_to_inmode, path_to_java, len_of_site, tmp_dir, "test"):
            true_scores.append(true_score)
        for false_score in calculate_scores_inmode_bootsrap(path_to_inmode, path_to_java, len_of_site, tmp_dir, "shuffled"):
            false_scores.append(false_score)
        shutil.rmtree(tmp_dir)
    table = creat_table_bootstrap(true_scores, false_scores)
    return(table)


def bootstrap_for_inmode(path, out, size_of_random_sample, path_to_inmode,
    order=2, path_to_java='java'):
    tmp_dir = os.getcwd() + '/tmp'
    sites = read_seqs(path)
    table = bootstrap_inmode(sites, path_to_inmode, path_to_java, tmp_dir, order, size_of_random_sample)
    write_table_bootstrap(out, table)
    return(0)

