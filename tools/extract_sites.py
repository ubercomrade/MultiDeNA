def read_sites(path):
    container = []
    with open(path) as file:
        for line in file:
            container.append(line.strip().split('\t')[-1])
    return(container)


def write_sites(sites, out_path):
    with open(out_path, 'w') as file:
        for line in sites:
            file.write(line + '\n')
    file.close()


def extract_sites(peaks_path, out_path):
    sites = read_sites(peaks_path)
    write_sites(sites, out_path)
    return(0)