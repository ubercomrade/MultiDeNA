import os
import sys
import subprocess
import argparse
import glob
import itertools
from operator import itemgetter
from shutil import copyfile
from tools.get_threshold_for_bamm import get_threshold_for_bamm
from tools.get_threshold_for_pwm import get_threshold_for_pwm
from tools.get_threshold_for_inmode import get_threshold_for_inmode
from tools.scan_by_pwm import scan_by_pwm
from tools.scan_by_bamm import scan_by_bamm
from tools.get_top_peaks import write_top_peaks
from tools.make_optimized_pwm import make_optimized_pwm
from tools.parse_chipmunk_results import parse_chipmunk_results
from tools.parse_inmode_results import parse_inmode_results
from tools.sites_intersection import sites_intersection
from tools.combine_results import combine_results

def prepare_data(path_to_genome, bed_path, bed, fasta, train_sample_size, test_sample_size):

    ########################
    #     GET TOP PEAKS    #
    ########################

    if not os.path.isfile(bed + '/' + 'train_sample.bed'):
        #Get top training_sample_size bed peaks
        print('Get top {0} bed peaks'.format(train_sample_size))
        bed_out = bed + '/'
        write_top_peaks(bed_path, bed_out, 4, 'train_sample', train_sample_size)      
    else:
        print('{0} already exists'.format('train_sample.bed'))

    if not os.path.isfile(bed + '/' + 'test_sample.bed'):
        #Get top testing_sample_size bed peaks
        print('Get top {0} bed peaks'.format(test_sample_size))
        bed_out = bed + '/'
        write_top_peaks(bed_path, bed_out, 4, 'test_sample', test_sample_size)
    else:
        print('{0} already exists'.format('test_sample.bed'))

    ########################
    #     BED TO FASTA     #
    ########################

    if not os.path.isfile(fasta + '/' + 'train_sample.fa'):
        #Bed peaks to fasta
        print('Get fasta from bed: {}'.format('train_sample.bed'))
        bed_to_fasta(path_to_genome,
            bed + '/train_sample.bed',
            fasta + '/train_sample.fa')
    else:
        print('{0} already exists'.format('train_sample.fa'))

    if not os.path.isfile(fasta + '/' + 'test_sample.fa'):
        print('Get fasta from bed: {}'.format('test_sample.bed'))
        bed_to_fasta(path_to_genome,
            bed + '/test_sample.bed',
            fasta + '/test_sample.fa')
    else:
        print('{0} already exists'.format('test_sample.fa'))
    return(0)


def get_inmode_model(models_path, fasta_path, path_to_java, path_to_inmode, motif_length, model_order):
    inmode_model_path = models_path + '/inmode_model'
    if not os.path.isdir(inmode_model_path):
        os.mkdir(inmode_model_path)
    if glob.glob(inmode_model_path + '/Learned_DeNovo*') == []:
        print('calculate INMODE model')
        args = [path_to_java, '-jar', path_to_inmode, 'denovo',
                'i={}'.format(fasta_path),
                'm={}'.format(motif_length),
               'mo={}'.format(model_order),
               'outdir={}'.format(inmode_model_path)]
        r = subprocess.call(args)
    else:
        print('INMODE model already exists')
    copyfile(glob.glob(inmode_model_path + '/Learned_DeNovo*/XML*')[0], inmode_model_path + '/inmode_model.xml')
    copyfile(glob.glob(inmode_model_path + '/Learned_DeNovo*/Binding_sites*')[0], inmode_model_path + '/inmode_sites.txt')
    return(0)


def run_chipmunk(path_to_java, path_to_chipmunk, fasta_path, path_out, motif_length_start, motif_length_end, cpu_count):
    args = [path_to_java, '-cp', path_to_chipmunk,
                   'ru.autosome.ChIPMunk', str(motif_length_start), str(motif_length_end), 'yes', '1.0',
                   's:{}'.format(fasta_path),
                  '100', '10', '1', str(cpu_count), 'random']
    p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
    out = p.communicate()
    with open(path_out, 'wb') as file:
        file.write(out[0])
    return(0)


def get_pwm_model(models_path, fasta_path,  path_to_java, path_to_chipmunk, motif_length_start, motif_length_end, cpu_count):
    chipmunk_model_path = models_path + '/pwm_model'
    if not os.path.isdir(chipmunk_model_path):
        os.mkdir(chipmunk_model_path)
    # FIND MODEL BY CHIPMUNK #
    if not os.path.isfile(chipmunk_model_path + '/initial_pwm_model.pwm'):
        print('Create PWM model by ChIPMunk')
        run_chipmunk(path_to_java, path_to_chipmunk,
        fasta_path,
        chipmunk_model_path + '/initial_pwm_model.txt',
        motif_length_start, motif_length_end, cpu_count)
    else:
        print('{0} already exists (initial model exists)'.format(chipmunk_model_path + '/initial_pwm_model.pwm'))
    # Parse results of chipmunk into files .meme, .pwm and .fasta (multi fasta) #
    parse_chipmunk_results(chipmunk_model_path + '/initial_pwm_model.txt',
        chipmunk_model_path, 'initial_pwm_model')
    # Get oPWM from chipmunk results. OUTPUT: .meme, .pwm and .fasta (multi fasta) #
    if not os.path.isfile(chipmunk_model_path + '/optimized_pwm_model.pwm'):
        make_optimized_pwm(chipmunk_model_path + '/initial_pwm_model.txt',
            fasta_path, chipmunk_model_path, 5000, 'optimized_pwm_model', cpu_count)
    else:
        print('{0} already exists'.format(chipmunk_model_path + '/optimized_pwm_model.pwm'))
    return(0)


def get_bamm_model(models_path, fasta_train, meme_model, model_order):
    #Get BaMM motif
    bamm_model_path = models_path + '/bamm_model'
    if not os.path.isfile(bamm_model_path + '/bamm_motif_1.ihbcp'):
        if not os.path.isdir(bamm_model_path):
            os.mkdir(bamm_model_path)
        args = ['BaMMmotif', bamm_model_path,
                fasta_train,
               '--PWMFile', meme_model,
                '--basename', 'bamm',
               '--EM',
               '--Order', str(model_order),
               '--order', str(model_order),
               '--scoreSeqset',
               '--saveLogOdds']
        r = subprocess.call(args)
    else:
        print('BAMM model already exists')
    return(0)


def calculate_thresholds_for_bamm(path_to_promoters, bamm_model_dir, thresholds_dir):
    if not os.path.isfile(thresholds_dir + '/bamm_model_thresholds.txt'):
        print('Calculate threshold for BAMM based on promoters and fpr')
        get_threshold_for_bamm(path_to_promoters,
            bamm_model_dir + '/bamm_motif_1.ihbcp',
            bamm_model_dir + '/bamm.hbcp',
            thresholds_dir + '/bamm_model_thresholds.txt')
    else:
        print('Thresholds for BAMM already calculated')
    return(0)


def calculate_thresholds_for_pwm(path_to_promoters, pwm_model_dir, thresholds_dir):
    if not os.path.isfile(thresholds_dir + '/pwm_model_thresholds.txt'):
        print('Calculate threshold for pwm based on promoters and fpr')
        get_threshold_for_pwm(path_to_promoters,
                pwm_model_dir + '/optimized_pwm_model.pwm',
                thresholds_dir + '/pwm_model_thresholds.txt')
    else:
        print('Thresholds for PWM already calculated')
    return(0)


def calculate_thresholds_for_inmode(path_to_promoters, inmode_model_dir, thresholds_dir, motif_length, path_to_inmode, path_to_java):
    if not os.path.isfile(thresholds_dir + '/inmode_model_thresholds.txt'):
        get_threshold_for_inmode(path_to_promoters,
                inmode_model_dir + '/inmode_model.xml',
                path_to_inmode,
                motif_length,
                thresholds_dir + '/inmode_model_thresholds.txt',
                path_to_java,
                thresholds_dir + '/tmp')
    else:
        print('Thresholds for INMODE already calculated')
    return(0)


def scan_peaks_by_pwm(fasta_test, model_path, scan, threshold_table_path, fpr):
    thr_pwm = get_threshold(threshold_table_path, fpr)
    pwm_scan_path = scan + '/pwm_{:.2e}.bed'.format(fpr)
    print('Scan peaks by PWM with FPR: {0} THR: {1}'.format(fpr, thr_pwm))
    scan_by_pwm(fasta_test, model_path, thr_pwm, pwm_scan_path)
    return(0)


def scan_peaks_by_bamm(fasta_test, model_path, bg_model_path, scan, threshold_table_path, fpr):
    thr_bamm = get_threshold(threshold_table_path, fpr)
    bamm_scan_path = scan + '/bamm_{:.2e}.bed'.format(fpr)
    print('Scan peaks by BAMM with FPR: {0} THR: {1}'.format(fpr, thr_bamm))
    scan_by_bamm(fasta_test, model_path, bg_model_path, thr_bamm, bamm_scan_path)
    return(0)


def scan_peaks_by_inmode(fasta_test, model_path, scan, threshold_table_path, fpr, path_to_java, path_to_inmode, path_to_promoters):
    inmode_scan_dir = scan + '/tmp'
    inmode_scan_path = scan + '/inmode_{:.2e}.bed'.format(fpr)
    thr_inmode = get_threshold(threshold_table_path, fpr)
    print('Scan peaks by INMODE with FPR: {0} THR: {1}'.format(fpr, thr_inmode))
    args = [path_to_java, '-Xmx6G', '-Xms1024m', '-jar', path_to_inmode, 'scan',
            'i={}'.format(model_path),
            'id={}'.format(fasta_test),
            'b={}'.format('From file'),
            'd={}'.format(path_to_promoters),
           'f={}'.format(fpr),
           'outdir={}'.format(inmode_scan_dir)]
    r = subprocess.call(args)
    parse_inmode_results(fasta_test, glob.glob(inmode_scan_dir + '/*.BED')[0],
        inmode_scan_path, thr_inmode)
    os.system("rm -r {}".format(inmode_scan_dir))
    return(0)


def get_sitega_model():
    pass


def bootstrap_sitega():
    pass


def calculate_thresholds_for_sitega():
    pass


def scan_peaks_by_sitega():
    pass


def bed_to_fasta(path_to_fa, path_to_bed, out):
    args = ['bedtools', 'getfasta' , '-s', '-name+',
            '-fi', path_to_fa,
            '-bed', path_to_bed,
            '-fo', out]
    r = subprocess.call(args)
    pass


def get_threshold(path, fpr_for_thr):
    container = list()
    append = container.append
    
    with open(path, 'r') as file:
        file.readline()
        for line in file:
            append(tuple(map(float, line.strip().split())))
    file.close()
    
    container = sorted(container, key=itemgetter(1))
    last_score, last_fpr = container[0]
    for line in container:
        if line[1] > fpr_for_thr:
            break
        else:
            last_score, last_fpr = line
    return(last_score)


def run_tomtom(query, model, outdir):
    args = ['tomtom', query, model, '-oc', outdir]
    r = subprocess.call(args)
    pass


def get_motif_length(models):
    with open(models + '/pwm_model/optimized_pwm_model.fasta', 'r') as file:
        for i in file:
            if i.startswith('>'):
                continue
            else:
                motif_length = len(i.strip())
                break
    file.close()
    return(motif_length)


# def compare_by_pair(bed, first, second, tag, name1, name2, compare_sites, path_to_python_tools):
#     list(itertools.combinations(a, 2))
#     args = ['pypy3', path_to_python_tools + '/compare_scripts/compare_sites_2.py',
#                     '-p', bed,
#                     '-first', first,
#                     '-second', second,
#                     '-t', tag, '-fname', name1, '-sname', name2,
#                     '-o', compare_sites]
#     r = subprocess.call(args)
#     pass


# def compare_2(bed, first, second, tag, name1, name2, compare_sites, path_to_python_tools):
#     args = ['pypy3', path_to_python_tools + '/compare_scripts/compare_sites_2.py',
#                     '-p', bed,
#                     '-first', first,
#                     '-second', second,
#                     '-t', tag, '-fname', name1, '-sname', name2,
#                     '-o', compare_sites]
#     r = subprocess.call(args)
#     pass


# def compare_4(bed, first, second, third, fourth, tag, compare_sites, path_to_python_tools):
#     args = ['python3', path_to_python_tools + '/compare_scripts/compare_sites_4.py',
#                     '-p', bed,
#                     '-first', first,
#                     '-second', second,
#                     '-third', third,
#                     '-fourth', fourth,
#                     '-t', tag,
#                     '-o', compare_sites]
#     r = subprocess.call(args)
#     pass


# def montecarlo( scores1, scores2, thr1, thr2, length, results):
#         args = ['monteCarlo', '{}'.format(scores1), '{}'.format(scores2), '{}'.format(thr1), '{}'.format(thr2), '{}'.format(length), '{}'.format(results)]
#         r = subprocess.call(args)
#         pass


def pipeline(tools, bed_path, fpr, train_sample_size, test_sample_size,
                      path_to_out, path_to_java, path_to_inmode, path_to_chipmunk,
                      path_to_promoters, path_to_genome, path_to_hocomoco, cpu_count):

    main_out = path_to_out
    model_order = 2
    cpu_count = cpu_count
    motif_length_start = str(8)
    motif_length_end = str(12)
    if not os.path.isdir(main_out):
        os.mkdir(main_out)
    models = main_out + '/models'
    bootstrap = models + '/bootstrap'
    thresholds = models + '/thresholds'
    fasta = main_out + '/fasta'
    bed = main_out + '/bed'
    scan = main_out + '/scan'
    scan_best = main_out + '/scan-best'
    results = main_out + '/results'
    tomtom = main_out + '/tomtom'
    montecarlo = main_out + '/montecarlo'
    
    ########################
    #      CREATE DIRS     #
    ########################

    if not os.path.isdir(models):
        os.mkdir(models)
    if not os.path.isdir(bootstrap):
        os.mkdir(bootstrap)
    if not os.path.isdir(thresholds):
        os.mkdir(thresholds)
    if not os.path.isdir(fasta):
        os.mkdir(fasta)
    if not os.path.isdir(bed):
        os.mkdir(bed)
    if not os.path.isdir(scan):
        os.mkdir(scan)
    if not os.path.isdir(scan_best):
        os.mkdir(scan_best)
    if not os.path.isdir(results):
        os.mkdir(results)
    if not os.path.isdir(tomtom):
        os.mkdir(tomtom)
    # if not os.path.isdir(montecarlo):
    #     os.mkdir(montecarlo)

    # PREPARE BED AND FASTA FILES #
    prepare_data(path_to_genome, bed_path, bed, fasta, train_sample_size, test_sample_size)

    fasta_train = fasta + '/train_sample.fa'
    fasta_test = fasta + '/test_sample.fa'
    bed_test = bed + '/test_sample.bed'
    bed_train = bed + '/train_sample.bed'

    # CALCULATE CHIPMUNK MODEL #
    pwm_model = models + '/pwm_model/optimized_pwm_model.pwm'
    pwm_threshold_table = thresholds + '/pwm_model_thresholds.txt'
    get_pwm_model(models, fasta_train,
        path_to_java, path_to_chipmunk,
        motif_length_start, motif_length_end,
        cpu_count)

    # BOOTSTRAP
    print('Run bootstrap for pwm model')
    bootstrap_for_pwm(models + '/pwm_model/optimized_pwm_model.fasta',
        bootstrap + '/pwm_model.tsv', 10000)
    check = check_bootstrap(bootstrap + '/pwm_model.tsv')
    if check < 0.0001:
        # THRESHOLD
        calculate_thresholds_for_pwm(path_to_promoters, models + '/pwm_model', thresholds)
        check = check_threshold_table(thresholds + '/pwm_model_thresholds.txt')
        if check < fpr:
            # SCAN
            scan_peaks_by_pwm(fasta_test, pwm_model, scan, pwm_threshold_table, fpr)
            scan_best_by_pwm(scan_best + '/pwm.scores.txt',
                 pwm_model,
                 fasta_test)
            extract_sites(scan + '/pwm_{:.2e}.bed'.format(fpr), tomtom + '/pwm.sites.txt')
            write_model(tomtom + '/pwm.sites.txt', tomtom, 'pwm')
            #run_tomtom(path_to_hocomoco, tomtom + '/pwm.meme', tomtom + '/pwm')
        else:
            tools.remove('pwm')
            print('PWM model has poor table with thresholds')
            print('GO to next model')
    else:
        tools.remove('pwm')
        print('PWM model is very weak (TPR = 0.5 -> FPR = {})'.format(check))
        print('GO to next model')


    # CALCULATE INMODE MODEL WITH EM ALG #
    motif_length = get_motif_length(models)
    inmode_model = models + '/inmode_model/inmode_model.xml'
    inmode_threshold_table = thresholds + '/inmode_model_thresholds.txt'

    get_inmode_model(models, fasta_train, path_to_java,
        path_to_inmode, motif_length, model_order)
    # BOOTSTRAP
    print('Run bootstrap for inmode model')
    bootstrap_for_inmode(models + '/inmode_model/inmode_sites.txt',
        bootstrap + "/inmode_model.tsv",
        10000,
        path_to_inmode, model_order, path_to_java)
    check = check_bootstrap(bootstrap + '/inmode_model.tsv')
    if check < 0.0001:
        # THRESHOLDS
        calculate_thresholds_for_inmode(path_to_promoters, models + '/inmode_model',
            thresholds, motif_length,
            path_to_inmode, path_to_java)
        check = check_threshold_table(thresholds + '/inmode_model_thresholds.txt')
        if check < fpr:
            # SCAN
            scan_peaks_by_inmode(fasta_test, inmode_model, scan, inmode_threshold_table,
            fpr, path_to_java, path_to_inmode, path_to_promoters)
            scan_best_by_inmode(scan_best + '/inmode.scores.txt',
                inmode_model,
                fasta_test,
                path_to_inmode, path_to_java)
            extract_sites(scan + '/inmode_{:.2e}.bed'.format(fpr), tomtom + '/inmode.sites.txt')
            write_model(tomtom + '/inmode.sites.txt', tomtom, 'inmode')
            #run_tomtom(path_to_hocomoco, tomtom + '/inmode.meme', tomtom + '/inmode')
        else:
            tools.remove('inmode')
            print('INMODE model has poor table with thresholds')
            print('GO to next model')
    else:
        tools.remove('inmode')
        print('INMODE model is very weak (TPR = 0.5 -> FPR = {} > 0.0001)'.format(check))
        print('GO to next model')


    # CALCULATE BAMM MODEL WITH EM ALG #
    meme_model = models + '/pwm_model/optimized_pwm_model.meme'
    bamm_threshold_table = thresholds + '/bamm_model_thresholds.txt'
    bamm_model = models + '/bamm_model/bamm_motif_1.ihbcp'
    bg_bamm_model = models + '/bamm_model/bamm.hbcp'

    get_bamm_model(models, fasta_train, meme_model, model_order)
    # BOOTSTRAP
    bootstrap_for_bamm(models + '/bamm_model/bamm_motif_1.occurrence',
        bootstrap + "/bamm_model.tsv", 10000, model_order)
    check = check_bootstrap(bootstrap + '/bamm_model.tsv')
    if check < 0.0001:
        # THRESHOLDS
        calculate_thresholds_for_bamm(path_to_promoters, models + '/bamm_model', thresholds)
        check = check_threshold_table(thresholds + '/bamm_model_thresholds.txt')
        if check < fpr:
            # SCAN
            scan_peaks_by_bamm(fasta_test, bamm_model, bg_bamm_model, scan, bamm_threshold_table, fpr)
            scan_best_by_bamm(scan_best + '/bamm.scores.txt',
                bamm_model,
                bg_bamm_model,
                fasta_test)
            extract_sites(scan + '/bamm_{:.2e}.bed'.format(fpr), tomtom + '/bamm.sites.txt')
            write_model(tomtom + '/bamm.sites.txt', tomtom, 'bamm')
            #run_tomtom(path_to_hocomoco, tomtom + '/bamm.meme', tomtom + '/bamm')
        else:
            tools.remove('bamm')
            print('BAMM model has poor table with thresholds')
            print('GO to next model')
    else:
        tools.remove('bamm')
        print('BAMM model is very weak (TPR = 0.5 -> FPR = {} > 0.0001)'.format(check))
        print('GO to next model')


    # CALCULATE SITEGA MODEL #
    sitega_model = models + '/sitega_model/sitega.txt'
    inmode_threshold_table = thresholds + '/sitega_model_thresholds.txt'
    get_sitega_model()
    # BOOTSTRAP
    bootstrap_sitega()
    # THRESHOLDS
    calculate_thresholds_for_sitega()
    scan_peaks_by_sitega()

    # COMPARE SITES #
    tools.remove('inmode')
    print('COMPARE SITES')
    pair_tools = list(itertools.combinations(tools, 2))
    for tool1, tool2 in pair_tools:
        tag = 'compare'
        scan1 = scan + '/{0}_{1:.2e}.bed'.format(tool1, fpr)
        scan2 = scan + '/{0}_{1:.2e}.bed'.format(tool2, fpr)
        sites_intersection(bed_test, scan1, scan2, tag, tool1, tool2, results)

    # COMBINE SCAN
    list_bed_path = [scan + '/{0}_{1:.2e}.bed'.format(i, fpr) for i in tools]
    list_path_fpr_table = [thresholds + '{}_model_thresholds.txt'.format(i) for i in tools]
    combine_results(fasta_test, list_bed_path, list_path_fpr_table, tools, results + '\combined_scan.pro')

    # CALCULATE SUMMARY

    # MONTECARLO #
    # for tool1, tool2 in pair_tools:
    #     thr1 = str(get_threshold(thresholds + '/{}_model_thresholds.txt'.format(tool1), fpr))
    #     thr2 = str(get_threshold(thresholds + '/{}_model_thresholds.txt'.format(tool2), fpr))
    #     scores1 = scan_best + '/{}.scores.txt'.format(tool1)
    #     scores2 = scan_best + '/{}.scores.txt'.format(tool2)
    #     name1, name2 = tool1, tool2
    #     results_montecarlo = montecarlo + '/montecarlo.results.{0}.{1}.{2:.2e}.txt'.format(tool1, tool2, fpr)
    #     montecarlo( scores1, scores2, thr1, thr2, length, results_montecarlo)



def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('bed', action='store', help='path to BED file')
    parser.add_argument('promoters', action='store', help='path to promoters fasta file')
    parser.add_argument('genome', action='store', help='path to genome fasta file')
    parser.add_argument('output', action='store', help='output dir')
    parser.add_argument('models', action='store', choices=['pwm', 'bamm', 'inmode', 'sitega'], metavar='N', nargs='+',
         help='list of models to use (pwm, bamm, inmode, sitega)')
    parser.add_argument('-t', '--train', action='store', type=int, dest='train_size',
                        required=False, default=500, help='size of training sample, by default size is equal to 500')
    parser.add_argument('-f', '--FPR', action='store', type=float, dest='fpr',
                        required=False, default=1.9*10**(-4), help='FPR, def=1.9*10^(-4)')
    parser.add_argument('-T', '--test', action='store', type=int, dest='test_size',
                        required=False, default=4000, help='size of testing sample, by default size is equal to 4000')
    parser.add_argument('-I', '--inmode', action='store', dest='inmode',
                        required=True, help='path to inmode')
    parser.add_argument('-J', '--java', action='store', dest='java',
                    required=False, default="java", help='path to Java')
    parser.add_argument('-c', '--chipmunk', action='store', dest='chipmunk',
                        required=True, help='path to chipmunk')
    parser.add_argument('-C', '--processes', action='store', type=int, dest='cpu_count',
                        required=False, default=4, help='Number of processes to use, default: 2')
    parser.add_argument('-H', '--hocomoco', action='store', dest='path_to_hocomoco',
                        required=False, help='path to HOCOMOCO database in meme format for TOMTOM')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():
    args = parse_args()
    bed_path = args.bed
    path_to_out = args.output
    train_sample_size = args.train_size
    test_sample_size = args.test_size
    fpr = args.fpr
    tools = args.models

    path_to_java = args.java
    path_to_chipmunk = args.chipmunk
    path_to_inmode = args.inmode
    path_to_promoters = args.promoters
    path_to_genome = args.genome
    path_to_hocomoco = args.path_to_hocomoco
    cpu_count = args.cpu_count

    pipeline(tools, bed_path, fpr, train_sample_size, test_sample_size,
                          path_to_out, path_to_java, path_to_inmode, path_to_chipmunk,
                          path_to_promoters, path_to_genome, path_to_hocomoco, cpu_count)

if __name__ == '__main__':
    main()
