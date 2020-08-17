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
    print(len(scores))
    if len(scores) > 3:
        sorted_peaks = sorted(peaks, key=lambda i: float(i[col]), reverse=True)
        for index, line in enumerate(sorted_peaks):
            line[3] = 'peaks_' + str(index)
        results = sorted_peaks[:amount]
    else:
        results = random.choices(peaks, k=amount)
        for index, line in enumerate(results):
            line[3] = 'peaks_' + str(index)
        results = sorted(results, key=lambda i: int(i[3].split('_')[1]), reverse=False)
    return(results)


def get_legths(data):
    l = list()
    for line in data:
        l.append(int(line[2]) - int(line[1]))
    return(l)


def split_peaks(peaks):
    peaks1 = []  # for train
    peaks2 = []  # for test
    for i in range(len(peaks)):
        if i % 2 == 0:
            peaks1.append(peaks[i])
        else:
            peaks2.append(peaks[i])
    return(peaks1, peaks2)


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
    peaks = get_top_peaks(peaks, amount, col)
    l = get_legths(peaks)
    write_length(output + '/' + tag + '.length.txt', l)
    with open(output + '/' + tag + '.bed', 'w') as file:
        for i in peaks:
            i[-1] = str(i[-1])
            file.write('\t'.join(i) + '\n')
    file.close()
    return(0)