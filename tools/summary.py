import itertools


def read_profile(path):
    container = {}
    with open(path) as file:
        for line in file:
            if line.startswith('>'):
                key = line.strip().split('\t')[-1]
                container[key] = []
            else:
                start, end, logfpr, strand, site, model, cluster = line.strip().split()
                container[key].append({
                    'start': int(start),
                    'end': int(end),
                    'log10fpr': float(logfpr),
                    'strand': strand,
                    'site': site,
                    'model': model,
                    'cluster': cluster
                })
    return(container)


def get_tools_combinations(profile):
    tools = []
    for key in profile.keys():
        for line in profile[key]:
            tool = line['model']
            if not tool in tools:
                tools.append(tool)
            else: continue
    container = []
    for i in range(1, len(tools) + 1):
        container += list(itertools.combinations(tools, i))
    combinations = []
    for i in container:
        combinations.append(set(i))
    return(combinations)


def peak_classifier(profile, combinations):
    statistics = dict()
    for i in range(len(combinations)):
        statistics[i] = 0
    for k in profile.keys():
        models_in_peak = set()
        for site in profile[k]:
            models_in_peak.add(site['model'])
        if models_in_peak != set():
            statistics[combinations.index(models_in_peak)] += 1
    return(statistics)


def write_peaks_classification(profile_path, write_path):
    profile = read_profile(profile_path)
    number_of_peaks = len([k for k in profile.keys()])
    combinations = get_tools_combinations(profile)
    statistics = peak_classifier(profile, combinations)
    header = []
    for i in combinations:
        header.append('x'.join(list(i)))
    header.append('NO_SITES')
    with open(write_path, 'w') as file:
        file.write('\t'.join(header) + '\n')
        file.write('\t'.join(map(str, statistics.values())) + '\t{}\n'.format(number_of_peaks))
    return(0)