import math
import os
import re


def read_fasta(path):
    fasta = list()
    with open(path, 'r') as file:
        for line in file:
            #print(line)
            if line.startswith('>'):
                line = line[1:].strip().split(':')
                record = dict()
                record['name'] = int(line[0].split('_')[1])
                record['chr'] = line[2]
                coordinates_strand = line[3]
                start, end = re.findall(r'\d*-\d*', coordinates_strand)[0].split('-')
                record['start'] = start
                record['end'] = end

                strand = re.findall(r'\(.\)', coordinates_strand[:-3])
                if not strand == []:
                    record['strand'] = strand[0].strip('()')
                else:
                    record['strand'] = '+'
            else:
                record['seq'] = line.strip().upper()
                fasta.append(record)
    file.close()
    return(fasta)


def complement(seq):
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').replace('N', 'n').upper()[::-1])


def read_inmode_bed(path):
    container = []
    with open(path) as file:
        for line in file:
            n, start, end, strand, score = line.split()
            container.append({
                'chr': '.',
                'start': int(start),
                'end': int(end),
                'id': int(n),
                'score': math.log(float(score), 10),
                'strand': strand,
                'site': '.'
            })
    return(container)


def write_bed(data, threshold, path):
    with open(path, 'w') as file:
        for line in data:
            if line['score'] >= threshold:
                line = list(line.values())
                line = [str(i) for i in line]
                file.write('\t'.join(line) + '\n')
            else:
                continue
    return(0)


def parse_inmode_results(input_fasta, input_bed, output, threshold):
    if os.path.getsize(input_bed) > 1:
        fasta = read_fasta(input_fasta)
        bed = read_inmode_bed(input_bed)
        for index, line in enumerate(bed):
            record = fasta[int(line['id'])]
            if line['strand'] == '-':
                line['site'] = complement(record['seq'][line['start']:line['end']])
            else:
                line['site'] = record['seq'][line['start']:line['end']]
            line['chr'] = record['chr']
            line['id'] = 'peaks_' + str(record['name'])
            line['start'] = int(line['start']) + int(record['start'])
            line['end'] = int(line['end']) + int(record['start'])
        write_bed(bed, threshold, output)
    else:
        with open(output, 'w') as file:
            print("INMODE model didn't fined sites")
            file.write('')
        file.close()
    return(0)
