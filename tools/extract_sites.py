import pandas as pd
import numpy as np


def read_bed_like_file(path):
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'name', 'score', 'strand', 'site'])
    return(df)


def write_sites(peaks, out_path):
    sites = list(peaks['site'][np.logical_not(pd.isna(peaks['site']))])
    with open(out_path, 'w') as file:
        for line in sites:
            file.write(line + '\n')
    file.close()


def extract_sites(peaks_path, out_path):
    peaks = read_bed_like_file(peaks_path)
    write_sites(peaks, out_path)
    return(0)

