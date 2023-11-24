#!/usr/bin/env python3

import os
import os.path
import sys
import subprocess
import argparse
import glob
import itertools
import shutil
import fnmatch
import pkg_resources
import pandas as pd
from operator import itemgetter
from shutil import copyfile
from multidena.tools.creat_optimized_pwm_model_chipmunk import de_novo_with_oprimization_pwm_chipmunk
from multidena.tools.creat_optimized_pwm_model_streme import de_novo_with_oprimization_pwm_streme
from multidena.tools.creat_optimized_dipwm_model import de_novo_with_oprimization_dipwm
from multidena.tools.creat_optimized_bamm_model import de_novo_with_oprimization_bamm
from multidena.tools.creat_optimized_inmode_model import de_novo_with_oprimization_inmode
from multidena.tools.creat_optimized_sitega_model import de_novo_with_oprimization_sitega
from multidena.tools.get_threshold_for_bamm import get_threshold_for_bamm
from multidena.tools.get_threshold_for_pwm import get_threshold_for_pwm
from multidena.tools.get_threshold_for_dipwm import get_threshold_for_dipwm
from multidena.tools.get_threshold_for_inmode import get_threshold_for_inmode
from multidena.tools.scan_by_pwm import scan_by_pwm
from multidena.tools.scan_by_dipwm import scan_by_dipwm
from multidena.tools.scan_by_bamm import scan_by_bamm
from multidena.tools.get_top_peaks import write_top_peaks
from multidena.tools.parse_chipmunk_results import parse_chipmunk_results
from multidena.tools.parse_inmode_results import parse_inmode_results
from multidena.tools.sites_intersection import sites_intersection
from multidena.tools.combine_results import combine_results_pro_format, combine_results_bed_format
from multidena.tools.summary import write_peaks_classification
from multidena.tools.scan_best_by_pwm import scan_best_by_pwm
from multidena.tools.scan_best_by_dipwm import scan_best_by_dipwm
from multidena.tools.scan_best_by_bamm import scan_best_by_bamm
from multidena.tools.scan_best_by_inmode import scan_best_by_inmode
from multidena.tools.extract_sites import extract_sites
from multidena.tools.write_model import write_model
from multidena.tools.clear_from_n import clear_from_n
from multidena.tools.parse_sitega_results import parse_sitega_results
from multidena.lib.common import check_threshold_table, check_bootstrap, gene_associated_with_motifs



# try:
#     from tools.creat_optimized_strum_model import de_novo_with_oprimization_strum
#     from tools.get_threshold_for_strum import get_threshold_for_strum
#     from tools.scan_by_strum import scan_by_strum
#     from tools.scan_best_by_strum import scan_best_by_strum
# except ModuleNotFoundError:
#     print('StruM module not found. You can`t use StruM model')
#     pass


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

        print('Get top {0} bed peaks'.format(test_sample_size))
        bed_out = bed + '/'
        write_top_peaks(bed_path, bed_out, 4, 'test_sample', test_sample_size)

        # if not os.path.isfile(bed + '/' + 'test_sample.bed'):
        #     #Get top testing_sample_size bed peaks
        #     print('Get top {0} bed peaks'.format(test_sample_size))
        #     bed_out = bed + '/'
        #     write_top_peaks(bed_path, bed_out, 4, 'test_sample', test_sample_size)
        # else:
        #     print('{0} already exists'.format('test_sample.bed'))

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


    print('Get fasta from bed: {}'.format('test_sample.bed'))
    bed_to_fasta(path_to_genome,
        bed + '/test_sample.bed',
        fasta + '/test_sample.fa')

    # if not os.path.isfile(fasta + '/' + 'test_sample.fa'):
    #     print('Get fasta from bed: {}'.format('test_sample.bed'))
    #     bed_to_fasta(path_to_genome,
    #         bed + '/test_sample.bed',
    #         fasta + '/test_sample.fa')
    # else:
    #     print('{0} already exists'.format('test_sample.fa'))
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
        print('Thresholds for diPWM already calculated')
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


# def calculate_thresholds_for_strum(path_to_promoters, strum_model_dir, thresholds_dir):
#     if not os.path.isfile(thresholds_dir + '/strum_model_thresholds.txt'):
#         print('Calculate threshold for StruM based on promoters and fpr')
#         get_threshold_for_strum(path_to_promoters,
#                 strum_model_dir + '/strum_model.pickle',
#                 thresholds_dir + '/strum_model_thresholds.txt')
#     else:
#         print('Thresholds for StruM already calculated')
#     return(0)


def scan_peaks_by_pwm(fasta_test, model_path, scan, threshold_table_path, fpr, tag):
    thr_pwm = get_threshold(threshold_table_path, fpr)
    pwm_scan_path = scan + '/pwm_{0}_{1:.2e}.bed'.format(tag, fpr)
    print('Scan peaks ({2}) by PWM with FPR: {0} THR: {1}'.format(fpr, thr_pwm, tag))
    scan_by_pwm(fasta_test, model_path, thr_pwm, pwm_scan_path)
    return(0)


def scan_peaks_by_dipwm(fasta_test, model_path, scan, threshold_table_path, fpr, tag):
    thr_pwm = get_threshold(threshold_table_path, fpr)
    dipwm_scan_path = scan + '/dipwm_{0}_{1:.2e}.bed'.format(tag, fpr)
    print('Scan peaks ({2}) by diPWM with FPR: {0} THR: {1}'.format(fpr, thr_pwm, tag))
    scan_by_dipwm(fasta_test, model_path, thr_pwm, dipwm_scan_path)
    return(0)


def scan_peaks_by_bamm(fasta_test, model_path, bg_model_path, scan, threshold_table_path, fpr, tag):
    thr_bamm = get_threshold(threshold_table_path, fpr)
    bamm_scan_path = scan + '/bamm_{0}_{1:.2e}.bed'.format(tag, fpr)
    print('Scan peaks ({2}) by BAMM with FPR: {0} THR: {1}'.format(fpr, thr_bamm, tag))
    scan_by_bamm(fasta_test, model_path, bg_model_path, thr_bamm, bamm_scan_path)
    return(0)


def scan_peaks_by_inmode(fasta_test, model_path, scan, threshold_table_path, fpr,
    path_to_java, path_to_inmode, path_to_promoters, tag):
    inmode_scan_dir = scan + '/tmp'
    inmode_scan_path = scan + '/inmode_{0}_{1:.2e}.bed'.format(tag, fpr)
    thr_inmode = get_threshold(threshold_table_path, fpr)
    print('Scan peaks ({2}) by INMODE with FPR: {0} THR: {1}'.format(fpr, thr_inmode, tag))
    args = [path_to_java, '-Xmx6G', '-Xms1024m', '-jar', path_to_inmode, 'scan',
            'i={}'.format(model_path),
            'id={}'.format(fasta_test),
            'b={}'.format('From file'),
            'd={}'.format(path_to_promoters),
           'f={}'.format(fpr),
           'outdir={}'.format(inmode_scan_dir)]
    #print(' '.join(args))
    r = subprocess.run(args, capture_output=True)
    out = r.stderr
    #print(out)
    parse_inmode_results(fasta_test, glob.glob(inmode_scan_dir + '/*.BED')[0],
        inmode_scan_path, thr_inmode)
    os.system("rm -r {}".format(inmode_scan_dir))
    return(0)


# def scan_peaks_by_strum(fasta_test, model_path, scan, threshold_table_path, fpr, tag):
#     thr_strum = get_threshold(threshold_table_path, fpr)
#     strum_scan_path = scan + '/strum_{0}_{1:.2e}.bed'.format(tag, fpr)
#     print('Scan peaks ({2}) by StruM with FPR: {0} THR: {1}'.format(fpr, thr_strum, tag))
#     scan_by_strum(fasta_test, model_path, thr_strum, strum_scan_path)
#     return(0)


def calculate_thresholds_for_sitega(path_to_promoters, sitega_model, thresholds_dir):
    dir_to_promoters = os.path.dirname(path_to_promoters)
    name_of_promoters = os.path.basename(path_to_promoters)
    if not os.path.isfile(thresholds_dir + '/sitega_model_thresholds.txt'):
        print('Calculate threshold for SiteGA based on promoters and fpr')
        args = ['sitega_thr_dist_mat',
            '{}/'.format(dir_to_promoters),
            '{}'.format(sitega_model),
            '{}'.format(name_of_promoters),
            '{}'.format(thresholds_dir + '/sitega_model_thresholds.txt'),
            '{}'.format(0.0005),
           '{}'.format(0.5),
           '{}'.format(0.0000000005)]
        r = subprocess.run(args, capture_output=True)
    else:
        print('Thresholds for SiteGA already calculated')
    return(0)


def scan_peaks_by_sitega(fasta_test, model_path, sitega_length, scan, threshold_table_path, fpr, scan_best_dir, tag):
    sitega_scan_tmp_dir = scan + '/sitega.tmp/'
    sitega_scan_path = scan + '/sitega_{0}_{1:.2e}.bed'.format(tag, fpr)
    thr_sitega = get_threshold(threshold_table_path, fpr)
    if not os.path.exists(sitega_scan_tmp_dir):
        os.mkdir(sitega_scan_tmp_dir)
    print('Scan peaks ({2}) by SiteGA with FPR: {0} THR: {1}'.format(fpr, thr_sitega, tag))
    args = ['andy1_mat',
        '{}'.format(fasta_test),
        '{}'.format(model_path),
        '{}'.format(threshold_table_path),
        '0',
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



# def get_motif_length(models):
#     with open(models + '/pwm_model/pwm_model.fasta', 'r') as file:
#         for i in file:
#             if i.startswith('>'):
#                 continue
#             else:
#                 motif_length = len(i.strip())
#                 break
#     file.close()
#     return(motif_length)


def run_annotation(list_of_scans, list_of_models, genome, output_dir):
    main_directory = pkg_resources.resource_filename('multidena', 'scripts')
    r_path = os.path.join(main_directory, 'annotation.R')
    args = [r_path,
        '--input_annotations', ';'.join(list_of_scans),
        '--models_names', ';'.join(list_of_models),
        '--genome', genome,
        '--output_dir', os.path.abspath(output_dir)]
    r = subprocess.run(args, capture_output=True)
    return(0)


def plot_motifs(dir_with_meme, output_dir):
    main_directory = pkg_resources.resource_filename('multidena', 'scripts')
    #r_path = os.path.join(main_directory,
    #    'scripts/motifCompare.R')
    r_path = os.path.join(main_directory,
        'motifCompare.R')
    args = [r_path,
        '--dir_with_motifs', dir_with_meme,
        '--dir_to_write', output_dir]
    r = subprocess.run(args, capture_output=True)
    return(0)


def get_length_sitega_model(path):
    with open(path) as file:
        file.readline()
        file.readline()
        length = int(file.readline().strip().split()[0])
    return(length)


def get_properties(path):
    with open(path) as file:
        length, order = file.readline().strip().split('\t')
    return int(length), int(order)


def get_peaks_with_sites(peaks_path, scan_path, write_path):
    peaks = pd.read_csv(peaks_path, sep='\t', header=None)
    if os.stat(scan_path).st_size != 0:
        scan = pd.read_csv(scan_path, sep='\t', header=None)
        scan_ids = set(scan[3])
        peaks_with_sites = peaks[[True if i in scan_ids else False for i in peaks[3]]]
        peaks_with_sites.to_csv(write_path, sep='\t', index=False, header=None)
    else:
        print(f'Empty file: {scan_path}')
        open(write_path, 'a').close()
    return 0


def pipeline(tools, bed_path, background_path, fpr, train_sample_size, test_sample_size,
                      path_to_out, path_to_java, path_to_inmode, path_to_chipmunk,
                      path_to_promoters, path_to_genome, organism, path_to_mdb):

    main_out = path_to_out
    model_order = 2
    cpu_count = 1
    pfpr = 0.001
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
    annotation = main_out + '/annotation'
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
    if not os.path.isdir(annotation):
        os.mkdir(annotation)

    # if not os.path.isdir(montecarlo):
    #     os.mkdir(montecarlo)

    # PREPARE BED AND FASTA FILES #
    prepare_data(path_to_genome, bed_path, bed, fasta, train_sample_size, test_sample_size)

    fasta_train = fasta + '/train_sample.fa'
    fasta_test = fasta + '/test_sample.fa'
    bed_test = bed + '/test_sample.bed'
    bed_train = bed + '/train_sample.bed'

    ### CALCULATE PWM MODEL ###
    if 'pwm-chipmunk' in tools or 'pwm-streme' in tools:
        pwm_model = models + '/pwm_model/pwm_model.pwm'
        pwm_threshold_table = thresholds + '/pwm_model_thresholds.txt'
        if not os.path.isfile(pwm_model):
            print('Training PWM model')
            if 'pwm-chipmunk' in tools:
                pwm_length = de_novo_with_oprimization_pwm_chipmunk(fasta_train, background_path, path_to_java, path_to_chipmunk,
                    models + '/pwm.tmp', models + '/pwm_model/',
                    output_auc + '/pwm', cpu_count, pfpr)
            else:
                pwm_length = de_novo_with_oprimization_pwm_streme(fasta_train,
                    background_path,
                    models + '/pwm.tmp/',
                    models + '/pwm_model/',
                    output_auc + '/pwm/',
                    pfpr)
        # THRESHOLD
        calculate_thresholds_for_pwm(path_to_promoters, models + '/pwm_model', thresholds)
        check = check_threshold_table(pwm_threshold_table)
        if check < fpr:
            # SCAN
            scan_peaks_by_pwm(fasta_train, pwm_model, scan, pwm_threshold_table, fpr, 'train')
            scan_peaks_by_pwm(fasta_test, pwm_model, scan, pwm_threshold_table, fpr, 'test')
            scan_best_by_pwm(scan_best + '/pwm.scores.txt',
                 pwm_model,
                 fasta_train)
            extract_sites(scan + '/pwm_train_{:.2e}.bed'.format(fpr), tomtom + '/pwm.sites.txt')
            write_model(tomtom + '/pwm.sites.txt', tomtom, 'pwm')
        else:
            print('WARNING! PWM model has poor table with thresholds')
            print('Best FPR for model is {}'.format(check))
            scan_peaks_by_pwm(fasta_train, pwm_model, scan, pwm_threshold_table, check, 'train')
            scan_best_by_pwm(scan_best + '/pwm.scores.txt',
                 pwm_model,
                 fasta_train)
            extract_sites(scan + '/pwm_train_{:.2e}.bed'.format(check), tomtom + '/pwm.sites.txt')
            write_model(tomtom + '/pwm.sites.txt', tomtom, 'pwm')
            os.remove(scan + '/pwm_train_{:.2e}.bed'.format(check))
            open(scan + '/pwm_train_{:.2e}.bed'.format(fpr), 'w').close()
            open(scan + '/pwm_test_{:.2e}.bed'.format(fpr), 'w').close()
    ### END PWM ###


    ### CALCULATE diPWM MODEL ###
    if 'dipwm' in tools:
        dipwm_model = models + '/dipwm_model/dipwm_model.pwm'
        dipwm_threshold_table = thresholds + '/dipwm_model_thresholds.txt'
        if not os.path.isfile(dipwm_model):
            print('Training diPWM model')
            dipwm_length = de_novo_with_oprimization_dipwm(fasta_train, background_path, path_to_java, path_to_chipmunk,
                models + '/dipwm.tmp', models + '/dipwm_model/',
                output_auc + '/dipwm', cpu_count, pfpr)
        # THRESHOLD
        calculate_thresholds_for_dipwm(path_to_promoters, models + '/dipwm_model', thresholds)
        check = check_threshold_table(dipwm_threshold_table)
        if check < fpr:
            # SCAN
            scan_peaks_by_dipwm(fasta_train, dipwm_model, scan, dipwm_threshold_table, fpr, 'train')
            scan_peaks_by_dipwm(fasta_test, dipwm_model, scan, dipwm_threshold_table, fpr, 'test')
            scan_best_by_dipwm(scan_best + '/dipwm.scores.txt',
                 dipwm_model,
                 fasta_train)
            extract_sites(scan + '/dipwm_train_{:.2e}.bed'.format(fpr), tomtom + '/dipwm.sites.txt')
            write_model(tomtom + '/dipwm.sites.txt', tomtom, 'dipwm')
        else:
            print('WARNING! diPWM model has poor table with thresholds')
            print('Best FPR for model is {}'.format(check))
            scan_peaks_by_dipwm(fasta_train, dipwm_model, scan, dipwm_threshold_table, check, 'train')
            scan_best_by_dipwm(scan_best + '/dipwm.scores.txt',
                 dipwm_model,
                 fasta_train)
            extract_sites(scan + '/dipwm_train_{:.2e}.bed'.format(check), tomtom + '/dipwm.sites.txt')
            write_model(tomtom + '/dipwm.sites.txt', tomtom, 'dipwm')
            os.remove(scan + '/dipwm_train_{:.2e}.bed'.format(check))
            open(scan + '/dipwm_train_{:.2e}.bed'.format(fpr), 'w').close()
            open(scan + '/dipwm_test_{:.2e}.bed'.format(fpr), 'w').close()
    ### END diPWM ###


    ### CALCULATE INMODE MODEL WITH EM ALG ###
    if 'inmode' in tools:
        #motif_length = get_motif_length(models)
        inmode_model_dir = models + '/inmode_model/'
        inmode_model = models + '/inmode_model/inmode_model.xml'
        inmode_threshold_table = thresholds + '/inmode_model_thresholds.txt'
        if not os.path.isfile(inmode_model):
            print('Training INMODE model')
            if not os.path.isdir(models + '/inmode_model/'):
                os.mkdir(models + '/inmode_model/')
            de_novo_with_oprimization_inmode(fasta_train, background_path, path_to_inmode,
                                            path_to_java, models + '/inmode.tmp',
                                            inmode_model_dir, output_auc + '/inmode',
                                            pfpr)
        # THRESHOLDS
        inmode_length, inmode_order = get_properties(f"{inmode_model_dir}/properties.txt")
        calculate_thresholds_for_inmode(path_to_promoters, models + '/inmode_model',
            thresholds, inmode_length,
            path_to_inmode, path_to_java)
        check = check_threshold_table(inmode_threshold_table)
        if check < fpr:
            # SCAN
            scan_peaks_by_inmode(fasta_train, inmode_model, scan, inmode_threshold_table,
            fpr, path_to_java, path_to_inmode, path_to_promoters, 'train')
            scan_peaks_by_inmode(fasta_test, inmode_model, scan, inmode_threshold_table,
            fpr, path_to_java, path_to_inmode, path_to_promoters, 'test')
            scan_best_by_inmode(scan_best + '/inmode.scores.txt',
                inmode_model,
                fasta_train,
                path_to_inmode, path_to_java, scan_best + '/inmode.tmp')
            extract_sites(scan + '/inmode_train_{:.2e}.bed'.format(fpr), tomtom + '/inmode.sites.txt')
            write_model(tomtom + '/inmode.sites.txt', tomtom, 'inmode')
        else:
            print('WARNING! INMODE model has poor table with thresholds')
            print('Best FPR for model is {}'.format(check))
            scan_peaks_by_inmode(fasta_train, inmode_model, scan, inmode_threshold_table,
            check, path_to_java, path_to_inmode, path_to_promoters, 'train')
            scan_best_by_inmode(scan_best + '/inmode.scores.txt',
                inmode_model,
                fasta_train,
                path_to_inmode, path_to_java, scan_best + '/inmode.tmp')
            extract_sites(scan + '/inmode_train_{:.2e}.bed'.format(check), tomtom + '/inmode.sites.txt')
            write_model(tomtom + '/inmode.sites.txt', tomtom, 'inmode')
            os.remove(scan + '/inmode_train_{:.2e}.bed'.format(check))
            open(scan + '/inmode_train_{:.2e}.bed'.format(fpr), 'w').close()
            open(scan + '/inmode_test_{:.2e}.bed'.format(fpr), 'w').close()
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
            bamm_length, bamm_order = de_novo_with_oprimization_bamm(fasta_train, background_path, output_auc + '/pwm',
                models + '/bamm.tmp', models + '/bamm_model', output_auc + '/bamm', pfpr)
        calculate_thresholds_for_bamm(path_to_promoters, models + '/bamm_model', thresholds)
        check = check_threshold_table(bamm_threshold_table)
        if check < fpr:
            # SCAN
            scan_peaks_by_bamm(fasta_train, bamm_model, bg_bamm_model, scan, bamm_threshold_table, fpr, 'train')
            scan_peaks_by_bamm(fasta_test, bamm_model, bg_bamm_model, scan, bamm_threshold_table, fpr, 'test')
            scan_best_by_bamm(scan_best + '/bamm.scores.txt',
                bamm_model,
                bg_bamm_model,
                fasta_train)
            extract_sites(scan + '/bamm_train_{:.2e}.bed'.format(fpr), tomtom + '/bamm.sites.txt')
            write_model(tomtom + '/bamm.sites.txt', tomtom, 'bamm')
        else:
            print('WARNING! BAMM model has poor table with thresholds')
            print('Best FPR for model is {}'.format(check))
            scan_peaks_by_bamm(fasta_train, bamm_model, bg_bamm_model, scan, bamm_threshold_table, check, 'train')
            scan_best_by_bamm(scan_best + '/bamm.scores.txt',
                bamm_model,
                bg_bamm_model,
                fasta_train)
            extract_sites(scan + '/bamm_train_{:.2e}.bed'.format(check), tomtom + '/bamm.sites.txt')
            write_model(tomtom + '/bamm.sites.txt', tomtom, 'bamm')
            os.remove(scan + '/bamm_train_{:.2e}.bed'.format(check))
            open(scan + '/bamm_train_{:.2e}.bed'.format(fpr), 'w').close()
            open(scan + '/bamm_test_{:.2e}.bed'.format(fpr), 'w').close()
    ### END BAMM ###


    ### CALCULATE SITEGA MODEL ###
    if 'sitega' in tools:
        sitega_model_dir =  models + '/sitega_model/'
        sitega_model_path = sitega_model_dir + 'sitega.mat'
        sitega_threshold_table = thresholds + '/sitega_model_thresholds.txt'
        if not os.path.isdir(sitega_model_dir):
            os.mkdir(sitega_model_dir)
        # PREPARE FASTA
        # clear_from_n(fasta_train, sitega_model_dir + '/train_sample_no_n.fa')
        # TRAIN SITEGA
        if not os.path.isfile(sitega_model_path):
            print('Training SiteGA model')
            de_novo_with_oprimization_sitega(
            fasta_train,
            background_path,
            models + '/sitega.tmp',
            models + '/sitega_model',
            output_auc + '/sitega',
            pfpr)
        sitega_length = get_length_sitega_model(sitega_model_path)
        calculate_thresholds_for_sitega(path_to_promoters, sitega_model_path, thresholds)
        check = check_threshold_table(sitega_threshold_table)
        if check < fpr:
            # SCAN
            scan_peaks_by_sitega(fasta_train, sitega_model_path, sitega_length,
                scan, sitega_threshold_table, fpr, scan_best, 'train')
            scan_peaks_by_sitega(fasta_test, sitega_model_path, sitega_length,
                scan, sitega_threshold_table, fpr, scan_best, 'test')
            extract_sites(scan + '/sitega_train_{:.2e}.bed'.format(fpr), tomtom + '/sitega.sites.txt')
            write_model(tomtom + '/sitega.sites.txt', tomtom, 'sitega')
        else:
            print('WARNING! SiteGA model has poor table with thresholds')
            print('Best FPR for model is {}'.format(check))
            scan_peaks_by_sitega(fasta_test, sitega_model_path, scan,
                sitega_threshold_table, check, scan_best , 'train')
            extract_sites(scan + '/sitega_train_{:.2e}.bed'.format(check), tomtom + '/sitega.sites.txt')
            write_model(tomtom + '/sitega.sites.txt', tomtom, 'sitega')
            os.remove(scan + '/sitega_train_{:.2e}.bed'.format(check))
            open(scan + '/sitega_train_{:.2e}.bed'.format(fpr), 'w').close()
            open(scan + '/sitega_test_{:.2e}.bed'.format(fpr), 'w').close()
    ### END SITEGA ###


    ### CALCULATE StruM MODEL ###
    # if 'strum' in tools:
    #     strum_model = models + '/strum_model/strum_model.pickle'
    #     strum_threshold_table = thresholds + '/strum_model_thresholds.txt'
    #     if not os.path.isfile(strum_model):
    #         print('Training StruM model')
    #         strum_length = de_novo_with_oprimization_strum(fasta_train, background_path,
    #             models + '/strum.tmp', models + '/strum_model/',
    #             output_auc + '/strum', cpu_count, pfpr)
    #     # THRESHOLD
    #     calculate_thresholds_for_strum(path_to_promoters, models + '/strum_model', thresholds)
    #     check = check_threshold_table(strum_threshold_table)
    #     if check < fpr:
    #         # SCAN
    #         scan_peaks_by_strum(fasta_train, strum_model, scan, strum_threshold_table, fpr, 'train')
    #         scan_peaks_by_strum(fasta_test, strum_model, scan, strum_threshold_table, fpr, 'test')
    #         scan_best_by_strum(scan_best + '/strum.scores.txt',
    #              strum_model,
    #              fasta_train)
    #         extract_sites(scan + '/strum_train_{:.2e}.bed'.format(fpr), tomtom + '/strum.sites.txt')
    #         write_model(tomtom + '/strum.sites.txt', tomtom, 'strum')
    #     else:
    #         print('WARNING! StruM model has poor table with thresholds')
    #         print('Best FPR for model is {}'.format(check))
    #         scan_peaks_by_strum(fasta_train, strum_model, scan, strum_threshold_table, check, 'train')
    #         scan_best_by_strum(scan_best + '/strum.scores.txt',
    #              strum_model,
    #              fasta_train)
    #         extract_sites(scan + '/strum_train_{:.2e}.bed'.format(check), tomtom + '/strum.sites.txt')
    #         write_model(tomtom + '/strum.sites.txt', tomtom, 'strum')
    #         os.remove(scan + '/strum_train_{:.2e}.bed'.format(check))
    #         open(scan + '/strum_train_{:.2e}.bed'.format(fpr), 'w').close()
    #         open(scan + '/strum_test_{:.2e}.bed'.format(fpr), 'w').close()
    ### END StruM ###


    # FIX NAME
    if 'pwm-chipmunk' in tools:
        tools[tools.index('pwm-chipmunk')] = 'pwm'
    else:
        tools[tools.index('pwm-streme')] = 'pwm'


    # COMPARE SITES
    print('Compare sites')
    #TRAIN
    pair_tools = list(itertools.combinations(tools, 2))
    for tool1, tool2 in pair_tools:
        tag = 'compare_train_{:.2e}'.format(fpr)
        scan1 = scan + '/{0}_{2}_{1:.2e}.bed'.format(tool1, fpr, 'train')
        scan2 = scan + '/{0}_{2}_{1:.2e}.bed'.format(tool2, fpr, 'train')
        sites_intersection(bed_train, scan1, scan2, tag, tool1, tool2, results)
    #TEST
    pair_tools = list(itertools.combinations(tools, 2))
    for tool1, tool2 in pair_tools:
        tag = 'compare_test_{:.2e}'.format(fpr)
        scan1 = scan + '/{0}_{2}_{1:.2e}.bed'.format(tool1, fpr, 'test')
        scan2 = scan + '/{0}_{2}_{1:.2e}.bed'.format(tool2, fpr, 'test')
        sites_intersection(bed_test, scan1, scan2, tag, tool1, tool2, results)


    # COMBINE SCAN
    for tag in ['train', 'test']:
        f = fasta + '/{}_sample.fa'.format(tag)
        list_bed_path = [scan + '/{0}_{1}_{2:.2e}.bed'.format(i, tag, fpr) for i in tools]
        list_path_fpr_table = [thresholds + '/{}_model_thresholds.txt'.format(i) for i in tools]
        combine_results_pro_format(f, list_bed_path, list_path_fpr_table,
            tools, results + '/combined_scan_{0}_{1:.2e}.pro'.format(tag, fpr))
        #combine_results_bed_format(f, list_bed_path, list_path_fpr_table, tools, results + '/combined_scan_{0}_{1:.2e}.bed'.format(tag, fpr))

        # CALCULATE SUMMARY
        write_peaks_classification(results + '/combined_scan_{0}_{1:.2e}.pro'.format(tag, fpr),
            tools, results + '/peaks_classification_{0}_{1:.2e}.tsv'.format(tag, fpr))


    # TOMTOM
    if os.path.exists(path_to_mdb):
        print('Runing TomTom')
        try:
            for t in tools:
                run_tomtom(path_to_mdb, tomtom + '/{}.meme'.format(t),
                    tomtom + '/{}.tomtom_results'.format(t))
        except:
            print('Tomtom failed. Check your PATH if you have already installed TomTom,\
                or install TomTom. If the problem rise again check your meme file (data base). \
                You can download motif database in meme format from http://meme-suite.org/doc/download.html.')


    # ANNOTATION AND GO
    # if organism in ['mm10', 'hg38', 'tair10']:
    #     name_converter = {
    #     'pwm': 'PWM',
    #     'dipwm': 'diPWM',
    #     'inmode': 'InMoDe',
    #     'bamm': 'BaMM',
    #     'sitega': 'SiteGA',
    #     'strum': 'StruM'
    #     }
    #     list_of_models = [name_converter[t] for t in tools]
    #     for tag in ['train', 'test']:
    #         output_dir = annotation + '/{}'.format(tag)
    #         if not os.path.isdir(output_dir):
    #             os.mkdir(output_dir)
    #         for tool in tools:
    #             write_path = f'{output_dir}/{tool}_{tag}_{fpr:.2e}.ann_genes.txt'
    #             path_scan = scan + '/{0}_{1}_{2:.2e}.bed'.format(tool, tag, fpr)
    #             path_ann = pkg_resources.resource_filename('multidena', f'promoters/{organism}.bed')
    #             r = gene_associated_with_motifs(path_scan, path_ann, write_path)
    #         list_of_ann = [f'{output_dir}/{tool}_{tag}_{fpr:.2e}.ann_genes.txt' for tool in tools]
    #         run_annotation(list_of_ann, list_of_models, organism, output_dir)

    if organism in ['mm10', 'hg38', 'tair10']:
        name_converter = {
        'pwm': 'PWM',
        'dipwm': 'diPWM',
        'inmode': 'InMoDe',
        'bamm': 'BaMM',
        'sitega': 'SiteGA',
        'strum': 'StruM'
        }
        list_of_models = [name_converter[t] for t in tools]
        for tag in ['train', 'test']:
            output_dir = annotation + '/{}'.format(tag)
            if not os.path.isdir(output_dir):
                os.mkdir(output_dir)
            for tool in tools:
                write_path = f'{output_dir}/{tool}_{tag}_{fpr:.2e}.ann_genes.txt'
                path_peaks = f'{output_dir}/{tool}_{tag}_{fpr:.2e}.peaks.txt'
                path_scan = scan + '/{0}_{1}_{2:.2e}.bed'.format(tool, tag, fpr)
                get_peaks_with_sites(bed + f'/{tag}_sample.bed', path_scan, path_peaks)
                path_ann = pkg_resources.resource_filename('multidena', f'promoters/{organism}.bed')
                r = gene_associated_with_motifs(path_peaks, path_ann, write_path)
            #All peaks intersection with promoters
            write_path = f'{output_dir}/{tag}_sample.ann_genes.txt'
            path_ann = pkg_resources.resource_filename('multidena', f'promoters/{organism}.bed')
            gene_associated_with_motifs(bed + f'/{tag}_sample.bed', path_ann, write_path)

            #GO
            list_of_ann = [f'{output_dir}/{tool}_{tag}_{fpr:.2e}.ann_genes.txt' for tool in tools]
            run_annotation(list_of_ann, list_of_models, organism, output_dir)

    # PLOT MOTIFS
    plot_motifs(tomtom, results)

    # FINISH
    print('Pipeline is finished!')
    tools = [t.upper() for t in tools]
    print('Results calculated for the next models:', *tools)
    return 0


def check_tools(tools, path_to_chipmunk, path_to_inmode):
    if 'inmode' in tools:
        if shutil.which("java") != None:
            out = subprocess.run(['java', '-version'], capture_output=True)
            fisrt, second = out.stderr.decode().split()[2].strip('\"').split('.')[:2]
            if f'{fisrt}.{second}' != '1.8' and 'inmode' in tools:
                print("You haven`t got Java 8 installed or Java 8 is not added to PATH. \
                You have to install Java 8 to use InMoDe \
                becouse InMoDe is crashed with using other version of Java ")
                print('inmode is removed from analisys.')
                tools.remove('inmode')

    if shutil.which("java") == None:
        if 'pwm-chipmunk' in tools or 'dipwm' in tools or 'inmode' in tools:
            print('There is no Java. You can`t use InMoDe, ChIPMunk(pwm) and diChIPMunk (dipwm) models. \
            As they require Java.')
            if 'pwm-chipmunk' in tools:
                print('pwm-chipmunk is replaced by pwm-streme.')
                tools[tools.index('pwm-chipmunk')] = 'pwm-streme'
            if 'dipwm' in tools:
                print('dipwm is removed from analisys.')
                tools.remove('dipwm')
            if 'inmode' in tools:
                print('inmode is removed from analisys.')
                tools.remove('inmode')

    if 'pwm-streme' in tools and shutil.which("streme") == None:
        print('You haven`t got Streme installed or Streme is not added to PATH. \
            You have to install Streme and add it to PATH to use pwm-streme. \
            You can visit - https://meme-suite.org/meme/')
        print('pwm-streme is removed from analisys.')
        tools.remove('pwm-streme')

    if 'bamm' in tools and shutil.which("BaMMmotif") == None:
        print('You haven`t got BaMMmotif installed or BaMMmotif is not added to PATH. \
            You have to install BaMMmotif and add it to PATH to use bamm. \
            You can visit - https://github.com/soedinglab/BaMMmotif2')
        print('bamm is removed from analisys.')
        tools.remove('bamm')

    if 'sitega' in tools and shutil.which("andy05.exe") == None:
        print('You haven`t got SiteGA (andy05.exe and andy0bsn5.exe) installed or SiteGA (andy05) is not added to PATH. \
            You have to install SiteGA (andy05.exe and andy0bsn5.exe) and add it to PATH to use sitega. \
            You can visit - https://github.com/parthian-sterlet/sitega')
        print('sitega is removed from analisys.')
        tools.remove('sitega')

    if 'dipwm' in tools and not os.path.exists(path_to_chipmunk):
        print('The path to ChIPMunk doesn`t exist.')
        print('dipwm is removed from analisys.')
        tools.remove('dipwm')

    if 'pwm-chipmunk' in tools and not os.path.exists(path_to_chipmunk):
        print('The path to ChIPMunk doesn`t exist.')
        print('pwm-chipmunk is removed from analisys.')
        tools.remove('pwm-chipmunk')

    if 'inmode' in tools and not os.path.exists(path_to_inmode):
        print('The path to InMoDe doesn`t exist.')
        print('inmode is removed from analisys.')
        tools.remove('inmode')

    if len(tools) < 2:
        sys.exit("Number of tools (models) is less than 2. Exit")

    return tools


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('bed', action='store', help='path to BED file')
    parser.add_argument('promoters', action='store', choices=['mm10', 'hg38', 'tair10', 'b73', 'dm6', 'ce235', 'r64'], metavar='N',
         help='Promoters of organism (hg38, mm10, tair10, b73)')
    parser.add_argument('genome', action='store', help='Path to genome fasta file')
    parser.add_argument('output', action='store', help='Output dir')
    parser.add_argument('-b', '--background', action='store', type=str, dest='background',
                        required=False, default='shuffled', help='Path to background. \
                        It is used for de novo and to estimate ROC and PRC. \
                        if it is not given background is generated by shuffling')
    parser.add_argument('models', action='store', choices=['pwm-chipmunk', 'pwm-streme', 'dipwm', 'bamm', 'inmode', 'sitega'],
                        metavar='M', nargs='+',
                        help='list of models to use (pwm-chipmunk, pwm-streme, dipwm, bamm, inmode, sitega)')
    parser.add_argument('-t', '--train', action='store', type=int, dest='train_size',
                        required=False, default=500, help='Number of peaks for training sample. The default value is 500')
    parser.add_argument('-f', '--FPR', action='store', type=float, dest='fpr',
                        required=False, default=1.9*10**(-4), help='FPR, def=1.9*10^(-4)')
    parser.add_argument('-T', '--test', action='store', type=int, dest='test_size',
                        required=False, default=-1, help='Number of peaks for testing sample. \
                        It could be any value starting from number of peaks used in traning sample. \
                        If parameter is -1 all peaks are used. The default value is -1.')
    parser.add_argument('-I', '--inmode', action='store', dest='inmode',
                        required=False, default='', help='Path to inmode (jar)')
    # parser.add_argument('-J', '--java', action='store', dest='java',
    #                 required=False, default="java", help='path to Java')
    parser.add_argument('-c', '--chipmunk', action='store', dest='chipmunk',
                        required=False, default='', help='Path to chipmunk (jar)')
    parser.add_argument('-m', '--motifdatabase', action='store', dest='path_to_mdb',
                        required=False, default='', help='Path to motif database in meme format for TOMTOM. \
                        You can get motif database from http://meme-suite.org/doc/download.html')

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
    background_path = args.background

    #path_to_java = args.java
    path_to_java = 'java'
    path_to_chipmunk = args.chipmunk
    path_to_inmode = args.inmode
    organism = args.promoters
    path_to_genome = args.genome
    path_to_mdb = args.path_to_mdb

    #Check tools
    check_tools(tools, path_to_chipmunk, path_to_inmode)

    if organism == 'mm10':
        path_to_promoters = pkg_resources.resource_filename('multidena', 'promoters/mm10.fasta')
    elif organism == 'hg38':
        path_to_promoters = pkg_resources.resource_filename('multidena', 'promoters/hg38.fasta')
    elif organism == 'tair10':
        path_to_promoters = pkg_resources.resource_filename('multidena', 'promoters/tair10.fasta')
    elif organism == 'b73':
        path_to_promoters = pkg_resources.resource_filename('multidena', 'promoters/b73_v5.fasta')
    elif organism == 'dm6':
        path_to_promoters = pkg_resources.resource_filename('multidena', 'promoters/dm6.fasta')
    elif organism == 'ce235':
        path_to_promoters = pkg_resources.resource_filename('multidena', 'promoters/ce235.fasta')
    elif organism == 'r64':
        path_to_promoters = pkg_resources.resource_filename('multidena', 'promoters/r64.fasta')

    pipeline(tools, bed_path, background_path, fpr, train_sample_size, test_sample_size,
                          path_to_out, path_to_java, path_to_inmode, path_to_chipmunk,
                          path_to_promoters, path_to_genome, organism, path_to_mdb)

if __name__ == '__main__':
    main()
