import os
import os.path
import sys
import subprocess
import argparse
import glob
import itertools
import shutil
from operator import itemgetter
from shutil import copyfile
from tools.creat_optimized_pwm_model import de_novo_with_oprimization_pwm
from tools.creat_optimized_bamm_model import de_novo_with_oprimization_bamm
from tools.creat_optimized_inmode_model import de_novo_with_oprimization_inmode
from tools.get_threshold_for_bamm import get_threshold_for_bamm
from tools.get_threshold_for_pwm import get_threshold_for_pwm
from tools.get_threshold_for_inmode import get_threshold_for_inmode
from tools.bootstrap_for_pwm import bootstrap_for_pwm
from tools.bootstrap_for_bamm import bootstrap_for_bamm
from tools.bootstrap_for_inmode import bootstrap_for_inmode
from tools.scan_by_pwm import scan_by_pwm
from tools.scan_by_bamm import scan_by_bamm
from tools.get_top_peaks import write_top_peaks
from tools.parse_chipmunk_results import parse_chipmunk_results
from tools.parse_inmode_results import parse_inmode_results
from tools.sites_intersection import sites_intersection
from tools.combine_results import combine_results
from tools.summary import write_peaks_classification
from tools.scan_best_by_pwm import scan_best_by_pwm
from tools.scan_best_by_bamm import scan_best_by_bamm
from tools.scan_best_by_inmode import scan_best_by_inmode
from tools.extract_sites import extract_sites
from tools.write_model import write_model
from tools.clear_from_n import clear_from_n
from tools.parse_sitega_results import parse_sitega_results
from lib.common import check_threshold_table, check_bootstrap

def prepare_data(path_to_genome, bed_path, bed, fasta, train_sample_size, test_sample_size):

    ########################
    #     GET TOP PEAKS    #
    ########################

    if train_sample_size == test_sample_size:
        if not os.path.isfile(bed + '/' + 'train_sample.bed'):
            print('Get top {0} bed peaks'.format(train_sample_size))
            bed_out = bed + '/'
            write_top_peaks(bed_path, bed_out, 4, 'train_sample', train_sample_size)
        else:
            print('{0} already exists'.format('train_sample.bed'))

        if not os.path.isfile(bed + '/' + 'test_sample.bed'):
            print('Get top {0} bed peaks'.format(test_sample_size))
            copyfile(bed + '/train_sample.bed', bed + '/test_sample.bed')
            copyfile(bed + '/train_sample.length.txt', bed + '/test_sample.length.txt')
        else:
            print('{0} already exists'.format('test_sample.bed')) 
    else:
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
        args = [path_to_java, '-jar', path_to_inmode, 'denovo',
                'i={}'.format(fasta_path),
                'm={}'.format(motif_length),
               'mo={}'.format(model_order),
               'outdir={}'.format(inmode_model_path)]
        r = subprocess.run(args, capture_output=True)
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
    p = subprocess.run(args, shell=False, capture_output=True)
    out = p.stdout
    with open(path_out, 'wb') as file:
        file.write(out)
    return(0)


def get_pwm_model(models_path, fasta_path,  path_to_java, path_to_chipmunk, motif_length_start, motif_length_end, cpu_count):
    chipmunk_model_path = models_path + '/pwm_model'
    if not os.path.isdir(chipmunk_model_path):
        os.mkdir(chipmunk_model_path)
    # FIND MODEL BY CHIPMUNK #
    if not os.path.isfile(chipmunk_model_path + '/initial_pwm_model.pwm'):
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
    if not os.path.isfile(chipmunk_model_path + '/pwm_model.pwm'):
        make_pwm(chipmunk_model_path + '/initial_pwm_model.txt',
            fasta_path, chipmunk_model_path, 5000, 'pwm_model', cpu_count)
    else:
        print('{0} already exists'.format(chipmunk_model_path + '/pwm_model.pwm'))
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
        r = subprocess.run(args, capture_output=True)
    else:
        print('BAMM model already exists')
    return(0)


def calculate_thresholds_for_bamm(path_to_promoters, bamm_model_dir, thresholds_dir):
    if not os.path.isfile(thresholds_dir + '/bamm_model_thresholds.txt'):
        print('Calculate threshold for BaMM based on promoters and fpr')
        get_threshold_for_bamm(path_to_promoters,
            bamm_model_dir + '/bamm_model.ihbcp',
            bamm_model_dir + '/bamm.hbcp',
            thresholds_dir + '/bamm_model_thresholds.txt')
    else:
        print('Thresholds for BAMM already calculated')
    return(0)


def calculate_thresholds_for_pwm(path_to_promoters, pwm_model_dir, thresholds_dir):
    if not os.path.isfile(thresholds_dir + '/pwm_model_thresholds.txt'):
        print('Calculate threshold for PWM based on promoters and fpr')
        get_threshold_for_pwm(path_to_promoters,
                pwm_model_dir + '/pwm_model.pwm',
                thresholds_dir + '/pwm_model_thresholds.txt')
    else:
        print('Thresholds for PWM already calculated')
    return(0)


def calculate_thresholds_for_inmode(path_to_promoters, inmode_model_dir, thresholds_dir, motif_length, path_to_inmode, path_to_java):
    if not os.path.isfile(thresholds_dir + '/inmode_model_thresholds.txt'):
        print('Calculate threshold for InMoDe based on promoters and fpr')
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
    r = subprocess.run(args, capture_output=True)
    parse_inmode_results(fasta_test, glob.glob(inmode_scan_dir + '/*.BED')[0],
        inmode_scan_path, thr_inmode)
    os.system("rm -r {}".format(inmode_scan_dir))
    return(0)


def get_sitega_model(models_dir, fasta_path):
    sitega_model_path = models_dir + '/sitega_model'
    if not os.path.isdir(sitega_model_path):
        os.mkdir(sitega_model_path)
    # FIND MODEL BY SITEGA
    clear_from_n(fasta_path, sitega_model_path + '/train_sample_no_n.fa')
    if not os.path.isfile(sitega_model_path + '/peaks.mnt'):
        args = ['monte0dg' ,'6', sitega_model_path + '/train_sample_no_n.fa', sitega_model_path + '/peaks.mnt']
        capture = subprocess.run(args, capture_output=True)
    if not os.path.isfile(sitega_model_path + '/train_sample_no_n.fa_mat'):
        args = ['andy02', sitega_model_path + '/peaks.mnt', '30', '10', '90', '10']
        capture = subprocess.run(args, capture_output=True)
    else:
        print('{0} already exists (initial model exists)'.format(sitega_model_path + '/train_sample.fa_mat'))
    pass


def bootstrap_sitega():
    pass


def calculate_thresholds_for_sitega(path_to_promoters, sitega_model_dir, thresholds_dir):
    if not os.path.isfile(thresholds_dir + '/sitega_model_thresholds.txt'):
        args = ['sitega_thr_dist_mat', \
                 sitega_model_dir + '/train_sample_no_n.fa_mat', \
                 path_to_promoters, thresholds_dir + '/sitega_model_thresholds.txt', \
                '0.0005', '0.997', '0.0000000005']
        capture = subprocess.run(args, capture_output=True)
    else:
        print('Thresholds for SITEGA already calculated')
    return(0)


def scan_peaks_by_sitega(fasta_test, sitega_model_dir, scan, threshold_table_path, fpr, scan_best_dir):
    thr_pwm = get_threshold(threshold_table_path, fpr)
    sitega_scan_path = scan + '/sitega_{:.2e}.bed'.format(fpr)
    print('Scan peaks by SITEGA with FPR: {0} THR: {1}'.format(fpr, thr_pwm))
    args = ['andy1_mat', fasta_test, sitega_model_dir + '/train_sample_no_n.fa.fa_mat', \
    threshold_table_path, str(fpr), sitega_model_dir + '/chipseq.pro']
    capture = subprocess.run(args, capture_output=True)
    parse_sitega_results(sitega_model_dir + '/chipseq.pro', sitega_scan_path)
    shutil.copyfile(fasta_test + '_bestscosg', scan_best_dir + '/sitega.scores.txt')
    return(0)


def bed_to_fasta(path_to_fa, path_to_bed, out):
    args = ['bedtools', 'getfasta' , '-s', '-name+',
            '-fi', path_to_fa,
            '-bed', path_to_bed,
            '-fo', out]
    r = subprocess.run(args, capture_output=True)
    pass


def get_threshold(path, fpr_for_thr):
    container = list()
    append = container.append
    with open(path, 'r') as file:
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


def run_tomtom(query, target, outdir):
    args = ['tomtom',
        query,
        target,
        '-thresh', '1.0',
        '-oc', outdir]
    r = subprocess.run(args, capture_output=True)
    return(1)   



def get_motif_length(models):
    with open(models + '/pwm_model/pwm_model.fasta', 'r') as file:
        for i in file:
            if i.startswith('>'):
                continue
            else:
                motif_length = len(i.strip())
                break
    file.close()
    return(motif_length)


def pipeline(tools, bed_path, fpr, train_sample_size, test_sample_size, bootstrap_flag,
                      path_to_out, path_to_java, path_to_inmode, path_to_chipmunk,
                      path_to_promoters, path_to_genome, path_to_mdb, cpu_count, tpr, pfpr):

    main_out = path_to_out
    model_order = 2
    cpu_count = cpu_count
    motif_length_start = str(8)
    motif_length_end = str(16)
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

    ### CALCULATE PWM MODEL ###
    if 'pwm' in tools:
        pwm_model = models + '/pwm_model/pwm_model.pwm'
        pwm_threshold_table = thresholds + '/pwm_model_thresholds.txt'
        if not os.path.isfile(pwm_model):
            print('Training PWM model')
            de_novo_with_oprimization_pwm(fasta_train, path_to_java, path_to_chipmunk, 
                './pwm.tmp', models + '/pwm_model/', cpu_count, tpr, pfpr)
        motif_length = get_motif_length(models)

        # BOOTSTRAP
        if bootstrap_flag and not os.path.isfile(bootstrap + '/pwm_model.tsv'):
            print('Run bootstrap for PWM model')
            bootstrap_for_pwm(fasta_train, bootstrap + '/pwm_model.tsv', \
                bootstrap + '/pwm_model_full.tsv', motif_length, \
                path_to_java, path_to_chipmunk, './pwm.tmp', cpu_count, counter=10000000)

        # THRESHOLD
        calculate_thresholds_for_pwm(path_to_promoters, models + '/pwm_model', thresholds)
        check = check_threshold_table(pwm_threshold_table)
        if check < fpr:
            # SCAN
            scan_peaks_by_pwm(fasta_test, pwm_model, scan, pwm_threshold_table, fpr)
            scan_best_by_pwm(scan_best + '/pwm.scores.txt',
                 pwm_model,
                 fasta_test)
            extract_sites(scan + '/pwm_{:.2e}.bed'.format(fpr), tomtom + '/pwm.sites.txt')
            write_model(tomtom + '/pwm.sites.txt', tomtom, 'pwm')
        else:
            tools.remove('pwm')
            print('PWM model has poor table with thresholds')
            print('GO to next model')
    ### END PWM ###


    ### CALCULATE INMODE MODEL WITH EM ALG ###
    if 'inmode' in tools:
        motif_length = get_motif_length(models)
        inmode_model = models + '/inmode_model/inmode_model.xml'
        inmode_threshold_table = thresholds + '/inmode_model_thresholds.txt'
        if not os.path.isfile(inmode_model):
            print('Training INMODE model')
            if not os.path.isdir(models + '/inmode_model/'):
                os.mkdir(models + '/inmode_model/')
            inmode_order = de_novo_with_oprimization_inmode(fasta_train, 
                motif_length, path_to_inmode, \
                path_to_java, './inmode.tmp', inmode_model, tpr, pfpr)
            with open(models + '/inmode_model/order.txt', 'w') as file:
                file.write(str(inmode_order))
            file.close()
        # BOOTSTRAP
        if bootstrap_flag and not os.path.isfile(bootstrap + '/inmode_model.tsv'):
            print('Run bootstrap for INMODE model')
            with open(models + '/inmode_model/order.txt') as file:
                inmode_order = int(file.readline().strip())
            file.close()
            bootstrap_for_inmode(fasta_train, bootstrap + '/inmode_model.tsv', \
                bootstrap + '/inmode_model_full.tsv', motif_length, \
                path_to_inmode, path_to_java, './inmode.tmp', counter=10000000, order=inmode_order)
        # THRESHOLDS
        calculate_thresholds_for_inmode(path_to_promoters, models + '/inmode_model',
            thresholds, motif_length,
            path_to_inmode, path_to_java)
        check = check_threshold_table(inmode_threshold_table)
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
        else:
            tools.remove('inmode')
            print('INMODE model has poor table with thresholds')
            print('GO to next model')
    ### END INMODE ###


    ### CALCULATE BAMM MODEL WITH EM ALG ###
    if 'bamm' in tools:
        meme_model = models + '/pwm_model/pwm_model.meme'
        bamm_threshold_table = thresholds + '/bamm_model_thresholds.txt'
        bamm_model = models + '/bamm_model/bamm_model.ihbcp'
        bg_bamm_model = models + '/bamm_model/bamm.hbcp'
        if not os.path.isfile(bamm_model):
            print('Training BAMM model')
            if not os.path.isdir(models + '/bamm_model/'):
                os.mkdir(models + '/bamm_model/')
            bamm_order = de_novo_with_oprimization_bamm(fasta_train, \
                motif_length, meme_model, './bamm.tmp', models + '/bamm_model',
                tpr, pfpr)
            with open(models + '/bamm_model/order.txt', 'w') as file:
                file.write(str(bamm_order))
            file.close()
        # BOOTSTRAP
        if bootstrap_flag and not os.path.isfile(bootstrap + '/bamm_model.tsv'):
            print('Run bootstrap for BAMM model')
            with open(models + '/bamm_model/order.txt') as file:
                bamm_order = int(file.readline().strip())
            file.close()
            bootstrap_for_bamm(fasta_train, bootstrap + '/bamm_model.tsv', \
                       bootstrap + '/bamm_model_full.tsv', motif_length, \
                       path_to_chipmunk, path_to_java, cpu_count, 
                       './bamm.tmp', counter = 10000000, order=bamm_order)
        calculate_thresholds_for_bamm(path_to_promoters, models + '/bamm_model', thresholds)
        check = check_threshold_table(bamm_threshold_table)
        if check < fpr:
            # SCAN
            scan_peaks_by_bamm(fasta_test, bamm_model, bg_bamm_model, scan, bamm_threshold_table, fpr)
            scan_best_by_bamm(scan_best + '/bamm.scores.txt',
                bamm_model,
                bg_bamm_model,
                fasta_test)
            extract_sites(scan + '/bamm_{:.2e}.bed'.format(fpr), tomtom + '/bamm.sites.txt')
            write_model(tomtom + '/bamm.sites.txt', tomtom, 'bamm')
        else:
            tools.remove('bamm')
            print('BAMM model has poor table with thresholds')
            print('GO to next model')
    ### END BAMM ###


    ### CALCULATE SITEGA MODEL ###
    if 'sitega' in tools:
        sitega_model = models + '/sitega_model'
        sitega_threshold_table = thresholds + '/sitega_model_thresholds.txt'
        if not os.path.isdir(sitega_model):
            os.mkdir(sitega_model)
        # PREPARE FASTA 
        clear_from_n(fasta_train, sitega_model + '/train_sample_no_n.fa')
        # TRAIN SITEGA
        print('Training SITEGA model')
        get_sitega_model(models, fasta_train)
        # BOOTSTRAP
        #print('Run bootstrap for SITEGA model')
        #if bootstrap_flag:
            #bootstrap_sitega()
        #else:
            #print("Bootstrap for SITEGA model already calculated -> PASS")
        #check = check_bootstrap(bootstrap + '/sitega_model.tsv')
            # THRESHOLDS
        calculate_thresholds_for_sitega(path_to_promoters, sitega_model, thresholds)
        check = check_threshold_table(sitega_threshold_table)
        if check < fpr:
            # SCAN
            scan_peaks_by_sitega(fasta_test, sitega_model, scan, sitega_threshold_table, fpr, scan_best)
            extract_sites(scan + '/sitega_{:.2e}.bed'.format(fpr), tomtom + '/sitega.sites.txt')
            write_model(tomtom + '/sitega.sites.txt', tomtom, 'sitega')
        else:
            tools.remove('sitega')
            print('SITEGA model has poor table with thresholds')
            print('GO to next model')
    ### END SITEGA ###


    # COMPARE SITES
    print('COMPARE SITES')
    pair_tools = list(itertools.combinations(tools, 2))
    for tool1, tool2 in pair_tools:
        tag = 'compare'
        scan1 = scan + '/{0}_{1:.2e}.bed'.format(tool1, fpr)
        scan2 = scan + '/{0}_{1:.2e}.bed'.format(tool2, fpr)
        sites_intersection(bed_test, scan1, scan2, tag, tool1, tool2, results)


    # COMBINE SCAN
    list_bed_path = [scan + '/{0}_{1:.2e}.bed'.format(i, fpr) for i in tools]
    list_path_fpr_table = [thresholds + '/{}_model_thresholds.txt'.format(i) for i in tools]
    combine_results(fasta_test, list_bed_path, list_path_fpr_table, tools, results + '/combined_scan.pro')


    # CALCULATE SUMMARY
    write_peaks_classification(results + '/combined_scan.pro', results + '/peaks_classification.tsv')


    # TOMTOM
    if path_to_mdb != None:   
        print('Runing TomTom')
        try:
            for t in tools:
                run_tomtom(path_to_mdb, tomtom + '/{}.meme'.format(t),
                    tomtom + '/{}.tomtom_results'.format(t))
        except:
            print('Tomtom failed. Check your PATH if you have already installed TomTom,\
                or install TomTom. If the problem rise again check your meme file (data base). \
                You can download motif database in meme format from http://meme-suite.org/doc/download.html.')
    print('Pipeline is finished!')
    tools = [t.upper() for t in tools]
    print('Results calculated for the next models:', *tools)
    

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('bed', action='store', help='path to BED file')
    parser.add_argument('promoters', action='store', choices=['mm10', 'hg38', 'tair10'], metavar='N',
         help='promoters of organism (hg38, mm10)')
    parser.add_argument('genome', action='store', help='path to genome fasta file')
    parser.add_argument('output', action='store', help='output dir')
    parser.add_argument('models', action='store', choices=['pwm', 'bamm', 'inmode', 'sitega'], metavar='N', nargs='+',
         help='list of models to use (pwm, bamm, inmode, sitega)')
    parser.add_argument('-t', '--train', action='store', type=int, dest='train_size',
                        required=False, default=2000, help='size of training sample, by default size is equal to 500')
    parser.add_argument('-f', '--FPR', action='store', type=float, dest='fpr',
                        required=False, default=1.9*10**(-4), help='FPR, def=1.9*10^(-4)')
    parser.add_argument('-T', '--test', action='store', type=int, dest='test_size',
                        required=False, default=2000, help='size of testing sample, by default size is equal to 4000')
    parser.add_argument('-b', '--bootstrap', action='store_true', dest='bootstrap',
                        required=False, help='Flag to calculate ROC for models, default is False')
    parser.add_argument('-I', '--inmode', action='store', dest='inmode',
                        required=True, help='path to inmode')
    parser.add_argument('-J', '--java', action='store', dest='java',
                    required=False, default="java", help='path to Java')
    parser.add_argument('-c', '--chipmunk', action='store', dest='chipmunk',
                        required=True, help='path to chipmunk')
    parser.add_argument('-C', '--processes', action='store', type=int, dest='cpu_count',
                        required=False, default=4, help='Number of processes to use, default: 2')
    parser.add_argument('-m', '--motifdatabase', action='store', dest='path_to_mdb',
                        required=False, default=None, help='path to motif database in meme format for TOMTOM. \
                        You can get motif database from http://meme-suite.org/doc/download.html')
    parser.add_argument('-tpr', '--lengthTPR', action='store', type=float, dest='tpr',
                        required=False, default=0.5, help='TECHNICAL, Calculate fpr at the tpr for choose optimal length of model')
    parser.add_argument('-pfpr', '--partionalFPR', action='store', dest='pfpr', type=float,
                        required=False, default=0.001, help='TECHNICAL, Threshold for calculating pAUC')

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
    bootstrap_flag = args.bootstrap

    path_to_java = args.java
    path_to_chipmunk = args.chipmunk
    path_to_inmode = args.inmode
    organism = args.promoters
    path_to_genome = args.genome
    path_to_mdb = args.path_to_mdb
    cpu_count = args.cpu_count

    tpr = args.tpr
    pfpr = args.pfpr

    this_dir, this_filename = os.path.split(__file__)
    if organism == 'mm10':
        path_to_promoters = os.path.join(this_dir, "promoters", "mm10.fasta")
    elif organism == 'hg38':
        path_to_promoters = os.path.join(this_dir, "promoters", "hg38.fasta")
    elif organism == 'tair10':
        path_to_promoters = os.path.join(this_dir, "promoters", "tair10.fasta")

    pipeline(tools, bed_path, fpr, train_sample_size, test_sample_size, bootstrap_flag,
                          path_to_out, path_to_java, path_to_inmode, path_to_chipmunk,
                          path_to_promoters, path_to_genome, path_to_mdb, cpu_count, tpr, pfpr)

if __name__ == '__main__':
    main()
