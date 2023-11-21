import os
import random


def read_file(path):
    peaks = []
    with open(path, 'r') as file:
        for line in file:
            if line.isspace():
                continue
            else:
                line = line.strip().split()
                peaks.append(line)
    file.close()
    return(peaks)


def clear_peaks(peaks):
    peaks = [i for i in peaks if len(i[0]) < 6 and i[0] != "chrMT"]
    return(peaks)

def get_top_peaks(peaks, amount, col):
    scores = list(set([i[col] for i in peaks]))
    if amount == -1:
        amount = len(peaks)
    elif len(peaks) < amount:
        print("Number of peaks ({0}) less than threshold amount ({1})".format(len(peaks),
            amount))
        amount = len(peaks)
    sorted_peaks = sorted(peaks, key=lambda i: float(i[col]), reverse=True)
    for index, line in enumerate(sorted_peaks):
        line[3] = 'peaks_' + str(index)
    results = sorted_peaks[:amount]
    return(results)


def get_random_peaks(peaks, amount):
    if amount > len(peaks):
        print(f'Number of peaks less than {amount}. Max number of peaks will be used instead')
        amount = len(peaks)
    if amount == -1:
        amount = len(peaks)
    results = random.sample(peaks, k=amount)
    key = len(peaks[0]) < 4
    for index, line in enumerate(results):
        if key:
            line.append('peaks_' + str(index))
        else:
            line[3] = 'peaks_' + str(index)
    results = sorted(results, key=lambda i: int(i[3].split('_')[1]), reverse=False)
    return(results)


def get_legths(data):
    l = list()
    for line in data:
        l.append(int(line[2]) - int(line[1]))
    return(l)


def write_length(path, data):
    with open(path, "w") as file:
        for line in data:
            file.write("{0}\n".format(line))
    file.close()
    return(0)


def write_top_peaks(path, output, col, tag, amount):
    if not os.path.isdir(output):
        os.mkdir(output)
    peaks = read_file(path)
    peaks = clear_peaks(peaks)
    if len(peaks[0]) >= 5:
        peaks = get_top_peaks(peaks, amount, col)
    else:
        peaks = get_random_peaks(peaks, amount)
    l = get_legths(peaks)
    write_length(output + '/' + tag + '.length.txt', l)
    with open(output + '/' + tag + '.bed', 'w') as file:
        for i in peaks:
            i[-1] = str(i[-1])
            file.write('\t'.join(i) + '\n')
    file.close()
    return(0)
