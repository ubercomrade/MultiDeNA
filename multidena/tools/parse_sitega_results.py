import argparse
import sys
import re
import csv


def parse_sitega(path, length):
    sitega = list()
    with open(path, 'r') as file:
        for line in file:
            #print(line)
            if line.startswith('>'):
                line = line[1:].strip().split(':')
                name = line[0]
                chromosome = line[2]
                coordinates_strand = line[3]
                start, end = re.findall(r'\d*-\d*', coordinates_strand)[0].split('-')
                start = int(start)
                end = int(end)
            else:
                record = dict()
                line = line.strip().split()
                site = line[3].upper()
                strand = line[2]
                score = float(line[1])
                left_pos = int(line[0])
                start_site = start + left_pos
                end_site = start + left_pos + length

                record['chr'] = chromosome
                record['start'] = start_site
                record['end'] = end_site
                record['name'] = name
                record['score'] = score
                record['strand'] = strand
                record['site'] = site
                sitega.append(record)
    file.close()
    return(sitega)


def write_bed(path_out, data):
    if len(data) == 0:
        open(path_out, 'w', newline='').close()
        print('WARNING! SITEGA didn`t find sites with current threshold')
    else:
        with open(path_out, 'w', newline='') as csvfile:
            fieldnames = data[0].keys()
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
            for row in data:
                writer.writerow(row)
    return(0)


def parse_sitega_results(path_in, path_out, length):
    sitega = parse_sitega(path_in, length)
    write_bed(path_out, sitega)
