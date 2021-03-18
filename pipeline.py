import os
import os.path
import sys
import subprocess
import argparse
import glob
import itertools
import shutil
import fnmatch
from operator import itemgetter
from shutil import copyfile
from tools.creat_optimized_pwm_model import de_novo_with_oprimization_pwm
from tools.creat_optimized_dipwm_model import de_novo_with_oprimization_dipwm
from tools.creat_optimized_bamm_model import de_novo_with_oprimization_bamm
from tools.creat_optimized_inmode_model import de_novo_with_oprimization_inmode
from tools.creat_optimized_strum_model import de_novo_with_oprimization_strum
from tools.get_threshold_for_bamm import get_threshold_for_bamm
from tools.get_threshold_for_pwm import get_threshold_for_pwm
from tools.get_threshold_for_dipwm import get_threshold_for_dipwm
from tools.get_threshold_for_inmode import get_threshold_for_inmode
from tools.get_threshold_for_strum import get_threshold_for_strum
# from tools.bootstrap_for_pwm import bootstrap_for_pwm
# from tools.bootstrap_for_dipwm import bootstrap_for_dipwm
# from tools.bootstrap_for_bamm import bootstrap_for_bamm
# from tools.bootstrap_for_inmode import bootstrap_for_inmode
# from tools.bootstrap_for_sitega import bootstrap_for_sitega
from tools.scan_by_pwm import scan_by_pwm
from tools.scan_by_dipwm import scan_by_dipwm
from tools.scan_by_bamm import scan_by_bamm
from tools.scan_by_strum import scan_by_strum
from tools.get_top_peaks import write_top_peaks
from tools.parse_chipmunk_results import parse_chipmunk_results
from tools.parse_inmode_results import parse_inmode_results
from tools.sites_intersection import sites_intersection
from tools.combine_results import combine_results_pro_format, combine_results_bed_format
from tools.summary import write_peaks_classification
from tools.scan_best_by_pwm import scan_best_by_pwm
from tools.scan_best_by_dipwm import scan_best_by_dipwm
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


def calculate_thresholds_for_dipwm(path_to_promoters, dipwm_model_dir, thresholds_dir):
    if not os.path.isfile(thresholds_dir + '/dipwm_model_thresholds.txt'):
        print('Calculate threshold for diPWM based on promoters and fpr')
        get_threshold_for_dipwm(path_to_promoters,
                dipwm_model_dir + '/dipwm_model.pwm',
                thresholds_dir + '/dipwm_model_thresholds.txt')
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


def calculate_thresholds_for_strum(path_to_promoters, pwm_model_dir, thresholds_dir):
    if not os.path.isfile(thresholds_dir + '/strum_model_thresholds.txt'):
        print('Calculate threshold for StruM based on promoters and fpr')
        get_threshold_for_strum(path_to_promoters,
                strum_model_dir + '/strum_model.pickle',
                thresholds_dir + '/strum_model_thresholds.txt')
    else:
        print('Thresholds for StruM already calculated')
    return(0)


def scan_peaks_by_pwm(fasta_test, model_path, scan, threshold_table_path, fpr):
    thr_pwm = get_threshold(threshold_table_path, fpr)
    pwm_scan_path = scan + '/pwm_{:.2e}.bed'.format(fpr)
    print('Scan peaks by PWM with FPR: {0} THR: {1}'.format(fpr, thr_pwm))
    scan_by_pwm(fasta_test, model_path, thr_pwm, pwm_scan_path)
    return(0)


def scan_peaks_by_dipwm(fasta_test, model_path, scan, threshold_table_path, fpr):
    thr_pwm = get_threshold(threshold_table_path, fpr)
    dipwm_scan_path = scan + '/dipwm_{:.2e}.bed'.format(fpr)
    print('Scan peaks by diPWM with FPR: {0} THR: {1}'.format(fpr, thr_pwm))
    scan_by_dipwm(fasta_test, model_path, thr_pwm, dipwm_scan_path)
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


def scan_peaks_by_strum(fasta_test, model_path, scan, threshold_table_path, fpr):
    thr_strum = get_threshold(threshold_table_path, fpr)
    strum_scan_path = scan + '/strum_{:.2e}.bed'.format(fpr)
    print('Scan peaks by StruM with FPR: {0} THR: {1}'.format(fpr, thr_strum))
    scan_by_strum(fasta_test, model_path, thr_strum, strum_scan_path)
    return(0)


def get_sitega_model(sitega_model_dir, sitega_length, fasta_path):
    if not os.path.isdir(sitega_model_dir):
        os.mkdir(sitega_model_dir)
    # FIND MODEL BY SITEGA
    clear_from_n(fasta_path, sitega_model_dir + '/train_sample_no_n.fa')
    if not os.path.isfile(sitega_model_dir + '/peaks.mnt'):
        args = ['monte0dg' ,'6', sitega_model_dir + '/train_sample_no_n.fa', sitega_model_dir + '/peaks.mnt']
        capture = subprocess.run(args, capture_output=True)
    if not os.path.isfile(sitega_model_dir + '/sitega.mat'):
        args = ['andy02', sitega_model_dir + '/peaks.mnt', str(sitega_length), '60', '90', '10']
        capture = subprocess.run(args, capture_output=True)
        for file in os.listdir(sitega_model_dir):
            if fnmatch.fnmatch(file, 'train_sample_no_n.fa_mat_*'):
                shutil.move('{0}/{1}'.format(sitega_model_dir, file), '{0}/sitega.mat'.format(sitega_model_dir))
        container = []
        with open('{0}/sitega.mat'.format(sitega_model_dir)) as file:
            for line in file:
                container.append(line)
        file.close()
        with open('{0}/sitega.mat'.format(sitega_model_dir), 'w') as file:
            file.write('sitega\n')
            for line in container[1:]:
                file.write(line)
    else:
        print('{0} already exists (initial model exists)'.format(sitega_model_dir + '/sitega.mat'))
    pass


def calculate_thresholds_for_sitega(path_to_promoters, sitega_model, thresholds_dir):
    if not os.path.isfile(thresholds_dir + '/sitega_model_thresholds.txt'):
        print('Calculate threshold for SiteGA based on promoters and fpr')
        args = ['sitega_thr_dist_mat',
            '{}'.format(sitega_model),
            '{}'.format(path_to_promoters),
            '{}'.format(thresholds_dir + '/sitega_model_thresholds.txt'),
            '{}'.format(0.0005),
           '{}'.format(0.995),
           '{}'.format(0.0000000005)]
        r = subprocess.run(args, capture_output=True)
    else:
        print('Thresholds for SiteGA already calculated')
    return(0)


def scan_peaks_by_sitega(fasta_test, model_path, sitega_length, scan, threshold_table_path, fpr, scan_best_dir):
    sitega_scan_tmp_dir = scan + '/sitega.tmp/'
    sitega_scan_path = scan + '/sitega_{:.2e}.bed'.format(fpr)
    thr_sitega = get_threshold(threshold_table_path, fpr)
    if not os.path.exists(sitega_scan_tmp_dir):
        os.mkdir(sitega_scan_tmp_dir)
    print('Scan peaks by SiteGA with FPR: {0} THR: {1}'.format(fpr, thr_sitega))
    args = ['andy1_mat',
        '{}'.format(fasta_test),
        '{}'.format(model_path),
        '{}'.format(threshold_table_path),
        '{}'.format(fpr),
       '{}'.format(sitega_scan_tmp_dir + '/sitega.pro')]
    r = subprocess.run(args, capture_output=True)
    parse_sitega_results(sitega_scan_tmp_dir + '/sitega.pro',
                         sitega_scan_path, sitega_length)
    shutil.copyfile(fasta_test + '_bestscosg', scan_best_dir + '/sitega.scores.txt')
    os.system("rm -r {}".format(sitega_scan_tmp_dir))
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


def pipeline(tools, bed_path, fpr, train_sample_size, test_sample_size,
                      path_to_out, path_to_java, path_to_inmode, path_to_chipmunk,
                      path_to_promoters, path_to_genome, path_to_mdb, cpu_count, pfpr):

    main_out = path_to_out
    model_order = 2
    cpu_count = cpu_count
    motif_length_start = str(8)
    motif_length_end = str(16)
    if not os.path.isdir(main_out):
        os.mkdir(main_out)
    models = main_out + '/models'
    #bootstrap = models + '/bootstrap'
    thresholds = models + '/thresholds'
    fasta = main_out + '/fasta'
    bed = main_out + '/bed'
    scan = main_out + '/scan'
    scan_best = main_out + '/scan-best'
    results = main_out + '/results'
    tomtom = main_out + '/tomtom'
    montecarlo = main_out + '/montecarlo'
    output_auc = main_out + '/auc'
    
    ########################
    #      CREATE DIRS     #
    ########################

    if not os.path.isdir(models):
        os.mkdir(models)
    # if not os.path.isdir(bootstrap):
    #     os.mkdir(bootstrap)
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
    if not os.path.isdir(output_auc):
        os.mkdir(output_auc)
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
            pwm_length = de_novo_with_oprimization_pwm(fasta_train, path_to_java, path_to_chipmunk, 
                models + '/pwm.tmp', models + '/pwm_model/', 
                output_auc + '/pwm', cpu_count, pfpr)
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
            print('WARNING! PWM model has poor table with thresholds')
            print('Best FPR for model is {}'.format(check))
            scan_peaks_by_pwm(fasta_test, pwm_model, scan, pwm_threshold_table, check)
            scan_best_by_pwm(scan_best + '/pwm.scores.txt',
                 pwm_model,
                 fasta_test)
            extract_sites(scan + '/pwm_{:.2e}.bed'.format(check), tomtom + '/pwm.sites.txt')
            write_model(tomtom + '/pwm.sites.txt', tomtom, 'pwm')
            os.remove(scan + '/pwm_{:.2e}.bed'.format(check))
            open(scan + '/pwm_{:.2e}.bed'.format(fpr), 'w').close()
    ### END PWM ###


    ### CALCULATE diPWM MODEL ###
    if 'dipwm' in tools:
        dipwm_model = models + '/dipwm_model/dipwm_model.pwm'
        dipwm_threshold_table = thresholds + '/dipwm_model_thresholds.txt'
        if not os.path.isfile(dipwm_model):
            print('Training diPWM model')
            dipwm_length = de_novo_with_oprimization_dipwm(fasta_train, path_to_java, path_to_chipmunk, 
                models + '/dipwm.tmp', models + '/dipwm_model/',
                output_auc + '/dipwm', cpu_count, pfpr)
        # THRESHOLD
        calculate_thresholds_for_dipwm(path_to_promoters, models + '/dipwm_model', thresholds)
        check = check_threshold_table(dipwm_threshold_table)
        if check < fpr:
            # SCAN
            scan_peaks_by_dipwm(fasta_test, dipwm_model, scan, dipwm_threshold_table, fpr)
            scan_best_by_dipwm(scan_best + '/dipwm.scores.txt',
                 dipwm_model,
                 fasta_test)
            extract_sites(scan + '/dipwm_{:.2e}.bed'.format(fpr), tomtom + '/dipwm.sites.txt')
            write_model(tomtom + '/dipwm.sites.txt', tomtom, 'dipwm')
        else:
            print('WARNING! diPWM model has poor table with thresholds')
            print('Best FPR for model is {}'.format(check))
            scan_peaks_by_dipwm(fasta_test, dipwm_model, scan, dipwm_threshold_table, check)
            scan_best_by_dipwm(scan_best + '/dipwm.scores.txt',
                 dipwm_model,
                 fasta_test)
            extract_sites(scan + '/dipwm_{:.2e}.bed'.format(check), tomtom + '/dipwm.sites.txt')
            write_model(tomtom + '/dipwm.sites.txt', tomtom, 'dipwm')
            os.remove(scan + '/dipwm_{:.2e}.bed'.format(check))
            open(scan + '/dipwm_{:.2e}.bed'.format(fpr), 'w').close()
    ### END diPWM ###


    ### CALCULATE INMODE MODEL WITH EM ALG ###
    if 'inmode' in tools:
        motif_length = get_motif_length(models)
        inmode_model_dir = models + '/inmode_model/'
        inmode_model = models + '/inmode_model/inmode_model.xml'
        inmode_threshold_table = thresholds + '/inmode_model_thresholds.txt'
        if not os.path.isfile(inmode_model):
            print('Training INMODE model')
            if not os.path.isdir(models + '/inmode_model/'):
                os.mkdir(models + '/inmode_model/')
            inmode_length, inmode_order = de_novo_with_oprimization_inmode(
                                            fasta_train, path_to_inmode, 
                                            path_to_java, models + '/inmode.tmp', 
                                            inmode_model_dir, output_auc + '/inmode', 
                                            pfpr)
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
                path_to_inmode, path_to_java, scan_best + '/inmode.tmp')
            extract_sites(scan + '/inmode_{:.2e}.bed'.format(fpr), tomtom + '/inmode.sites.txt')
            write_model(tomtom + '/inmode.sites.txt', tomtom, 'inmode')
        else:
            print('WARNING! INMODE model has poor table with thresholds')
            print('Best FPR for model is {}'.format(check))
            scan_peaks_by_inmode(fasta_test, inmode_model, scan, inmode_threshold_table,
            check, path_to_java, path_to_inmode, path_to_promoters)
            scan_best_by_inmode(scan_best + '/inmode.scores.txt',
                inmode_model,
                fasta_test,
                path_to_inmode, path_to_java, scan_best + '/inmode.tmp')
            extract_sites(scan + '/inmode_{:.2e}.bed'.format(check), tomtom + '/inmode.sites.txt')
            write_model(tomtom + '/inmode.sites.txt', tomtom, 'inmode')
            os.remove(scan + '/inmode_{:.2e}.bed'.format(check))
            open(scan + '/inmode_{:.2e}.bed'.format(fpr), 'w').close()
    ### END INMODE ###


    ### CALCULATE BAMM MODEL WITH EM ALG ###
    if 'bamm' in tools:
        meme_model = models + '/pwm_model/pwm_model.meme'
        bamm_threshold_table = thresholds + '/bamm_model_thresholds.txt'
        bamm_model_dir = models + '/bamm_model/'
        bamm_model = models + '/bamm_model/bamm_model.ihbcp'
        bg_bamm_model = models + '/bamm_model/bamm.hbcp'
        if not os.path.isfile(bamm_model):
            print('Training BAMM model')
            if not os.path.isdir(models + '/bamm_model/'):
                os.mkdir(models + '/bamm_model/')
            bamm_length, bamm_order = de_novo_with_oprimization_bamm(fasta_train, output_auc + '/pwm', 
                models + '/bamm.tmp', models + '/bamm_model', output_auc + '/bamm', pfpr)
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
            print('WARNING! BAMM model has poor table with thresholds')
            print('Best FPR for model is {}'.format(check))
            scan_peaks_by_bamm(fasta_test, bamm_model, bg_bamm_model, scan, bamm_threshold_table, check)
            scan_best_by_bamm(scan_best + '/bamm.scores.txt',
                bamm_model,
                bg_bamm_model,
                fasta_test)
            extract_sites(scan + '/bamm_{:.2e}.bed'.format(check), tomtom + '/bamm.sites.txt')
            write_model(tomtom + '/bamm.sites.txt', tomtom, 'bamm')
            os.remove(scan + '/bamm_{:.2e}.bed'.format(check))
            open(scan + '/bamm_{:.2e}.bed'.format(fpr), 'w').close()
    ### END BAMM ###


    ### CALCULATE SITEGA MODEL ###
    if 'sitega' in tools:
        sitega_length = 60
        sitega_model_dir =  models + '/sitega_model/'
        sitega_model_path = sitega_model_dir + 'sitega.mat'
        sitega_threshold_table = thresholds + '/sitega_model_thresholds.txt'
        if not os.path.isdir(sitega_model_dir):
            os.mkdir(sitega_model_dir)
        # PREPARE FASTA 
        clear_from_n(fasta_train, sitega_model_dir + '/train_sample_no_n.fa')
        # TRAIN SITEGA
        if not os.path.isfile(sitega_model_path):
            print('Training SiteGA model')
            get_sitega_model(sitega_model_dir, sitega_length, fasta_train)
        with open('{0}/sitega.mat'.format(sitega_model_dir)) as file:
            file.readline()
            lpd = int(file.readline().strip().split()[0])
            sitega_length = int(file.readline().strip().split()[0])
        file.close()
        # BOOTSTRAP
        if bootstrap_flag and not os.path.isfile(bootstrap + '/sitega_model.tsv'):
            print('Run bootstrap for SiteGA model')
            bootstrap_for_sitega(sitega_model_dir + '/train_sample_no_n.fa', 
                bootstrap + '/sitega_model.tsv',
                bootstrap + '/sitega_model_full.tsv',
                sitega_length, lpd, bootstrap + '/sitega.tmp/', counter=5000000)
        calculate_thresholds_for_sitega(path_to_promoters, sitega_model_path, thresholds)
        check = check_threshold_table(sitega_threshold_table)
        if check < fpr:
            # SCAN
            scan_peaks_by_sitega(fasta_test, sitega_model_path, sitega_length,
                scan, sitega_threshold_table, fpr, scan_best)
            extract_sites(scan + '/sitega_{:.2e}.bed'.format(fpr), tomtom + '/sitega.sites.txt')
            write_model(tomtom + '/sitega.sites.txt', tomtom, 'sitega')
        else:
            print('WARNING! SiteGA model has poor table with thresholds')
            print('Best FPR for model is {}'.format(check))
            scan_peaks_by_sitega(fasta_test, sitega_model_path, scan, 
                sitega_threshold_table, check, scan_best)
            extract_sites(scan + '/sitega_{:.2e}.bed'.format(check), tomtom + '/sitega.sites.txt')
            write_model(tomtom + '/sitega.sites.txt', tomtom, 'sitega')
            os.remove(scan + '/sitega_{:.2e}.bed'.format(check))
            open(scan + '/sitega_{:.2e}.bed'.format(fpr), 'w').close()
    ### END SITEGA ###


    ### CALCULATE StruM MODEL ###
    if 'strum' in tools:
        strum_model = models + '/strum_model/strum_model.pickle'
        strum_threshold_table = thresholds + '/strum_model_thresholds.txt'
        if not os.path.isfile(strum_model):
            print('Training StruM model')
            strum_length = de_novo_with_oprimization_strum(fasta_train, 
                models + '/strum.tmp', models + '/strum_model/', 
                output_auc + '/strum', cpu_count, pfpr)
        # THRESHOLD
        calculate_thresholds_for_strum(path_to_promoters, models + '/strum_model', thresholds)
        check = check_threshold_table(strum_threshold_table)
        if check < fpr:
            # SCAN
            scan_peaks_by_strum(fasta_test, strum_model, scan, strum_threshold_table, fpr)
            # scan_best_by_strum(scan_best + '/strum.scores.txt',
            #      strum_model,
            #      fasta_test)
            extract_sites(scan + '/strum_{:.2e}.bed'.format(fpr), tomtom + '/strum.sites.txt')
            write_model(tomtom + '/strum.sites.txt', tomtom, 'strum')
        else:
            print('WARNING! StruM model has poor table with thresholds')
            print('Best FPR for model is {}'.format(check))
            scan_peaks_by_strum(fasta_test, strum_model, scan, strum_threshold_table, check)
            # scan_best_by_strum(scan_best + '/strum.scores.txt',
            #      strum_model,
            #      fasta_test)
            extract_sites(scan + '/strum_{:.2e}.bed'.format(check), tomtom + '/strum.sites.txt')
            write_model(tomtom + '/strum.sites.txt', tomtom, 'pwm')
            os.remove(scan + '/strum_{:.2e}.bed'.format(check))
            open(scan + '/strum_{:.2e}.bed'.format(fpr), 'w').close()
    ### END StruM ###


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
    combine_results_pro_format(fasta_test, list_bed_path, list_path_fpr_table, tools, results + '/combined_scan.pro')
    #combine_results_bed_format(fasta_test, list_bed_path, list_path_fpr_table, tools, results + '/combined_scan.bed')

    # CALCULATE SUMMARY
    write_peaks_classification(results + '/combined_scan.pro', tools, results + '/peaks_classification.tsv')

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
    parser.add_argument('models', action='store', choices=['pwm', 'dipwm', 'bamm', 'inmode', 'sitega', 'strum'], metavar='N', nargs='+',
         help='list of models to use (pwm, dipwm, bamm, inmode, sitega, strum)')
    parser.add_argument('-t', '--train', action='store', type=int, dest='train_size',
                        required=False, default=2000, help='size of training sample, by default size is equal to 500')
    parser.add_argument('-f', '--FPR', action='store', type=float, dest='fpr',
                        required=False, default=1.9*10**(-4), help='FPR, def=1.9*10^(-4)')
    parser.add_argument('-T', '--test', action='store', type=int, dest='test_size',
                        required=False, default=2000, help='size of testing sample, by default size is equal to 4000')
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

    path_to_java = args.java
    path_to_chipmunk = args.chipmunk
    path_to_inmode = args.inmode
    organism = args.promoters
    path_to_genome = args.genome
    path_to_mdb = args.path_to_mdb
    cpu_count = args.cpu_count

    pfpr = args.pfpr

    this_dir, this_filename = os.path.split(__file__)
    if organism == 'mm10':
        path_to_promoters = os.path.join(this_dir, "promoters", "mm10.fasta")
    elif organism == 'hg38':
        path_to_promoters = os.path.join(this_dir, "promoters", "hg38.fasta")
    elif organism == 'tair10':
        path_to_promoters = os.path.join(this_dir, "promoters", "tair10.fasta")

    pipeline(tools, bed_path, fpr, train_sample_size, test_sample_size,
                          path_to_out, path_to_java, path_to_inmode, path_to_chipmunk,
                          path_to_promoters, path_to_genome, path_to_mdb, cpu_count, pfpr)

if __name__ == '__main__':
    main()
