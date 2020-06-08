import subprocess
import os
import shutil
import random
from lib.common import read_bamm, write_table_bootstrap
from lib.speedup import calculate_scores_bamm_bootstrap, creat_random_sample, creat_table_bootstrap


def read_log_odds_zoops(path):
    container = []
    append = container.append
    with open(path) as file:
        file.readline()
        for line in file:
            site = line.split()[4]
            if not 'N' in site:
                append(site)
            else:
                continue
    return(container)


def create_bamm_model(directory, order):
    fasta_path = directory + '/train.fasta'
    sites_path = directory + '/sites.txt'
    args = ['BaMMmotif', directory, fasta_path, '--bindingSiteFile', sites_path, '--EM', '--order', str(order), '--Order', str(order)]
    subprocess.call(args)
    bamm_path = directory + '/train_motif_1.ihbcp'
    bg_path = directory + '/train.hbcp'
    bamm, order = read_bamm(bamm_path, bg_path)
    return(bamm, order)


def create_train_fasta_and_sites(sites, directory):
    with open(directory + '/train.fasta', 'w') as file:
        for index, site in enumerate(sites):
            file.write('>{0}\n{1}\n'.format(index, site))
    with open(directory + '/sites.txt', 'w') as file:
        for site in sites:
            file.write('{}\n'.format(site))
    return(0)


def bootstrap_bamm(sites, size_of_random_sample, order, directory):
    true_scores = []
    false_scores = []
    number_of_sites = len(sites)
    length_of_site = len(sites[0])

    for i in range(10):
        train_sample = random.choices(sites, k=round(0.9 * number_of_sites))
        test_sample = [site for site in sites  if not site in train_sample]
        random_sample = creat_random_sample(test_sample, size_of_random_sample)
        create_train_fasta_and_sites(train_sample, directory)
        bamm, order = create_bamm_model(directory, order)
        for true_score in calculate_scores_bamm_bootstrap(test_sample, bamm, order, length_of_site):
            true_scores.append(true_score)
        for false_score in calculate_scores_bamm_bootstrap(random_sample, bamm, order, length_of_site):
            false_scores.append(false_score)
    table = creat_table_bootstrap(true_scores, false_scores)
    return(table)


def bootstrap_for_bamm(path, out, size_of_random_sample, order):
    tmp_dir = os.getcwd() + '/tmp'
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    sites = read_log_odds_zoops(path)
    table = bootstrap_bamm(sites, size_of_random_sample, order, tmp_dir)
    write_table_bootstrap(out, table)
    shutil.rmtree(tmp_dir)
    return(0)

