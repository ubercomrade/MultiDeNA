import random
from lib.common import read_seqs, sites_to_pwm, write_table_bootstrap
from lib.speedup import calculate_scores_pwm_bootstrap, creat_table_bootstrap, creat_random_sample


def bootstrap_pwm(sites, size_of_random_sample):
    true_scores = []
    false_scores = []
    number_of_sites = len(sites)
    len_of_site = len(sites[0])

    for i in range(10):
        train_sample = random.choices(sites, k=round(0.9 * number_of_sites))
        test_sample = [site for site in sites  if not site in train_sample]
        random_sample = creat_random_sample(test_sample, size_of_random_sample)
        pwm = sites_to_pwm(train_sample)
        for true_score in calculate_scores_pwm_bootstrap(test_sample, pwm):
            true_scores.append(true_score)
        for false_score in calculate_scores_pwm_bootstrap(random_sample, pwm):
            false_scores.append(false_score)
    table = creat_table_bootstrap(true_scores, false_scores)
    return(table)


def bootstrap_for_pwm(path, out, size_of_random_sample):
    sites = read_seqs(path)
    table = bootstrap_pwm(sites, size_of_random_sample)
    write_table_bootstrap(out, table)
    return(0)
