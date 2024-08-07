import os
import subprocess
import shutil
import math
from multidena.lib.common import read_seqs_with_complement, get_threshold


def calculate_scores_inmode_thresholds(path_to_inmode, path_to_model, path_to_fasta, path_to_java, tmp_dir):
    container = list()
    append = container.append
    args = [path_to_java, '-Xmx8096m', '-Xms1024m', '-jar', path_to_inmode, 'scan',
        'i={}'.format(path_to_model),
        'id={}'.format(path_to_fasta),
        'b={}'.format('From file'),
        'd={}'.format(path_to_fasta),
       'f={}'.format(0.005),
       'outdir={}'.format(tmp_dir)]
    r = subprocess.run(args, capture_output=True)
    with open(tmp_dir + "/Motif_hits_from_SequenceScan(0.005).BED") as file:
        for line in file:
            append(math.log10(float(line.strip().split()[4])))
    return(container)


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
