import subprocess
import glob
import math
import shutil
import os
from operator import itemgetter


def inmode_scan(path_to_inmode, path_java, input_data, input_model, tmp_dir,
                     fpr_for_thr=1):
    args = [path_java, '-Xmx6G', '-Xms1024m', '-jar', path_to_inmode, 'scan',
            'i={}'.format(input_model),
            'id={}'.format(input_data),
           'f={}'.format(fpr_for_thr),
           'outdir={}'.format(tmp_dir)]
    r = subprocess.run(args, capture_output=True)
    return()


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


def parse_inmode_results(inmode_bed, out):
    container = []
    table = read_inmode_bed(inmode_bed)
    table.sort(key=itemgetter(0, 4))
    last_index = 0
    for line in table:
        index = line[0]
        score = line[4]
        if last_index != index:
            container.append(last_score)
        last_score = score
        last_index = index
    container.append(score)
    with open(out, 'w') as file:
        for i in container:
            file.write('{}\n'.format(math.log(float(i), 10)))
    file.close()
    return(container)


def scan_best_by_inmode(out, path_to_model, fasta_path, path_to_inmode, path_to_java):
    tmp_dir = os.getcwd() + '/tmp'
    inmode_scan(path_to_inmode, path_to_java, fasta_path, path_to_model, tmp_dir)
    inmode_bed = glob.glob(tmp_dir + '/*.BED')[0]
    parse_inmode_results(inmode_bed, out)
    shutil.rmtree(tmp_dir)
    return(0)