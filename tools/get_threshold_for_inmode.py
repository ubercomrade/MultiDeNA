import os
import subprocess
import shutil
import math
from lib.common import read_seqs_with_complement


def calculate_scores_inmode_thresholds(path_to_inmode, path_to_model, path_to_fasta, path_to_java, tmp_dir):
    container = list()
    append = container.append
    args = [path_to_java, '-Xmx8096m', '-Xms1024m', '-jar', path_to_inmode, 'scan',
        'i={}'.format(path_to_model),
        'id={}'.format(path_to_fasta),
        'b={}'.format('From file'),
        'd={}'.format(path_to_fasta),
       'f={}'.format(0.001),
       'outdir={}'.format(tmp_dir)]
    r = subprocess.run(args, capture_output=True)
    with open(tmp_dir + "/Motif_hits_from_SequenceScan(0.001).BED") as file:
        for line in file:
            append(math.log10(float(line.strip().split()[4])))
    return(container)


def get_threshold(scores, number_of_sites, path_out):
    scores.sort(reverse=True) # big -> small
    with open(path_out, "w") as file:
        last_score = scores[0]
        for count, score in enumerate(scores[1:], 1):
            if score == last_score:
                continue
            elif count/number_of_sites > 0.0005:
                file.write("{0}\t{1}\n".format(last_score, count/number_of_sites))
                break
            elif score != last_score:
                file.write("{0}\t{1}\n".format(last_score, count/number_of_sites))
                last_score = score 
    file.close()
    return(0)


def get_threshold_for_inmode(path_to_fasta, path_to_model, path_to_inmode,
    length_of_site, path_out, path_to_java='java', tmp_dir='./tmp'):
    tmp_dir = os.getcwd() + '/tmp'
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    peaks = read_seqs_with_complement(path_to_fasta)
    number_of_sites = sum([len(range(len(peak) - length_of_site + 1)) for peak in peaks])
    scores = calculate_scores_inmode_thresholds(path_to_inmode, path_to_model, path_to_fasta, path_to_java, tmp_dir)
    get_threshold(scores, number_of_sites, path_out)
    shutil.rmtree(tmp_dir)
    return(0)