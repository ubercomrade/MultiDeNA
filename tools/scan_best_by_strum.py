import numpy as np
from strum import strum
from lib.common import read_strum

def read_fasta(path):
    fasta = list()
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                fasta.append(line.strip())
    file.close()
    return(fasta)


def write_list(path, data):
    with open(path, "w") as file:
        for line in data:
            file.write("{0}\n".format(line))
    file.close()
    pass


def scan_best_by_strum(results_path, strum_path, fasta_path):
    fasta = read_fasta(fasta_path)
    strum_model = read_strum(strum_path)
    results = []
    for seq in fasta:
        rseq = strum.rev_comp(seq)
        scores_1 = strum_model.score_seq(seq)
        scores_2 = strum_model.score_seq(rseq)
        i1 = np.argmax(scores_1)
        i2 = np.argmax(scores_2)
        if scores_1[i1] > scores_2[i2]:
            results.append(scores_1[i1])
        else:
            results.append(scores_2[i2])
    write_list(results_path, results)
    return(0)
