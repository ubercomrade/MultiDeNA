import subprocess
import glob
import math
import shutil


def inmode_scan(path_to_inmode, path_java, input_data, input_model, tmp_dir,
                     fpr_for_thr=1):
    args = [path_java, '-Xmx6G', '-Xms1024m', '-jar', path_to_inmode, 'scan',
            'i={}'.format(input_model),
            'id={}'.format(input_data),
           'f={}'.format(fpr_for_thr),
           'outdir={}'.format(tmp_dir)]
    r = subprocess.call(args)
    return()


def parse_inmode_results(inmode_bed, out):
    container = list()
    append = container.append
    df = pd.read_csv(inmode_bed, sep='\t', header=None)
    df = df.sort_values(by=[0,4], ascending=[True, True])
    last_index = 0
    for index, score in zip(df[0], df[4]):
        if last_index != index:
            append(last_score)
        last_score = score
        last_index = index
    append(score)
    with open(out, 'w') as file:
        for i in container:
            file.write('{}\n'.format(math.log(float(i), 10)))
    file.close()
    return(0)
    

def scan_best_by_inmode(out, path_to_model, fasta_path, path_to_inmode, path_to_java):
    tmp_dir = os.getcwd() + '/tmp'
    inmode_scan(path_to_inmode, path_to_java, fasta_path, path_to_model, tmp_dir)
    inmode_bed = glob.glob(tmp_dir + '/*.BED')[0]
    parse_inmode_results(inmode_bed, out)
    shutil.rmtree(tmp_dir)
    return(0)