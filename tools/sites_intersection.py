import csv
from operator import itemgetter


def read_bed(path):
    container = []
    with open(path) as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        for row in csvreader:
            container.append({
                'chr': row[0], 'start': int(row[1]), 'end': int(row[2])})
    container = sorted(container, key = itemgetter('chr', 'start'))
    return(container)


def read_scan(path):
    container = []
    with open(path) as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        for row in csvreader:
            container.append({
                'chr': row[0], 'start': int(row[1]), 'end': int(row[2]),
                'name': row[3], 'score':row[4], 'strand': row[5], 'site': row[6]})
    container = sorted(container, key = itemgetter('chr', 'start'))
    return(container)

    
def check_sites_intersection_in_peak(sites1, sites2):
    intersected_sites1 = []
    intersected_sites2 = []
    for i in sites1:
        for j in sites2:
            if j['start'] < i['end'] and j['end'] > i['start']:
                intersected_sites1.append(i)
                break
    not_intersected_sites1 = [i for i in sites1 if not i in intersected_sites1]
    for i in sites2:
        for j in sites1:
            if j['start'] < i['end'] and j['end'] > i['start']:
                intersected_sites2.append(i) 
                break
    not_intersected_sites2 = [i for i in sites2 if not i in intersected_sites2]
    return(intersected_sites1, not_intersected_sites1,
          intersected_sites2, not_intersected_sites2)
        

def get_indexes(peaks, sites):
    container = dict()
    chrs = list(set([i['chr'] for i in peaks]))
    chrs.sort()
    index = 0
    for chr_ in chrs:
        sub_peaks = [i for i in peaks if i['chr'] == chr_]
        sub_sites = [i for i in sites if i['chr'] == chr_]
        if len(sub_sites) == 0:
            index += len(sub_peaks)
            continue
        for p in sub_peaks:
            for ss in sub_sites:
                if p['start'] < ss['end'] and p['end'] > ss['start']:
                    if not index in container:
                        container[index] = [ss]
                    else:
                        container[index].append(ss)
            index += 1
    return(container)


def check_common_peaks(sites1, sites2, common_ids):
    container = []
    for id_ in common_ids:
        sites1_with_id = sites1[id_]
        sites2_with_id = sites2[id_]
        intersected_sites1, not_intersected_sites1, intersected_sites2, not_intersected_sites2 = check_sites_intersection_in_peak(sites1_with_id, sites2_with_id)
        container.append({
            'intersected_sites1': intersected_sites1,
            'not_intersected_sites1': not_intersected_sites1,
            'intersected_sites2': intersected_sites2,
            'not_intersected_sites2': not_intersected_sites2
        })
    return(container)


def get_not_intersect_sites(sites1, sites2):
    chrs = list(set([i['chr'] for i in sites1]))
    chrs.sort()
    not_intersected = []
    for chr_ in chrs:
        intersected = []
        sub_sites1 = [i for i in sites1 if i['chr'] == chr_]
        sub_sites2 = [i for i in sites2 if i['chr'] == chr_]
        if len(sub_sites2) == 0:
            continue
        for s1 in sub_sites1:
            for s2 in sub_sites2:
                if s1['start'] < s2['end'] and s1['end'] > s2['start']:
                    intersected.append(s1)
    
        not_intersected += [i for i in sub_sites1 if not i in intersected]
    return(not_intersected)


def write_table(path, data):
    with open(path, 'w') as csvfile:
        fieldnames = list(data.keys())
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerow(data)
    return(0)


def sites_intersection(peaks_path, first_model_path, second_model_path, tag, fname, sname, out_dir):
    peaks = read_bed(peaks_path)
    first_model_sites = read_scan(first_model_path)
    first_model_sites_dict = get_indexes(peaks, first_model_sites)
    ids_first_model = set(first_model_sites_dict.keys())

    second_model_sites = read_scan(second_model_path)
    second_model_sites_dict = get_indexes(peaks, second_model_sites)
    ids_second_model = set(second_model_sites_dict.keys())

    only_first_model_ids = ids_first_model - ids_second_model
    only_second_model_ids = ids_second_model - ids_first_model
    common_ids = ids_first_model & ids_second_model

    c = check_common_peaks(first_model_sites_dict, second_model_sites_dict, common_ids)

    statistics = {'{0}: against {1}'.format(fname, sname): len(only_first_model_ids),
                  '{0}: against {1}'.format(sname, fname): len(only_second_model_ids),
                  'overlapped:{},{}'.format(fname, sname): 0,
                  #'OVERLAPPED:++{0},{1}'.format(fname, sname): 0,
                  #'OVERLAPPED:++{0},{1}'.format(sname, fname): 0,
                  'not_overlapped:{},{}'.format(fname, sname): 0,
                  'no_sites:{},{}'.format(fname, sname): 0,
                  'number_of_peaks': len(peaks)}

    statistics['overlapped:{},{}'.format(fname, sname)] = len([i for i in c if len(i['intersected_sites1']) > 0])
    statistics['not_overlapped:{},{}'.format(fname, sname)] = len([i for i in c if len(i['intersected_sites1']) == 0])
    #statistics['OVERLAPPED:++{0},{1}'.format(fname, sname)] = len([i for i in c if len(i['intersected_sites1']) > 0 and len(i['not_intersected_sites1']) > 0])
    #statistics['OVERLAPPED:++{0},{1}'.format(sname, fname)] = len([i for i in c if len(i['intersected_sites2']) > 0 and len(i['not_intersected_sites2']) > 0])
    statistics['no_sites:{},{}'.format(fname, sname)] = len(peaks) - len(common_ids) - statistics['{0}: against {1}'.format(fname, sname)] - statistics['{0}: against {1}'.format(sname, fname)]
    write_table(out_dir + '/' + tag + '_{0}.{1}_counts.tsv'.format(fname, sname), statistics)
    return(0)
