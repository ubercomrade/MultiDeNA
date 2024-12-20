# MultiDeNA is pipeline for integration different models of transcription factor binding sites

## Introduction
Regulation of eukaryotic gene expression achieved through compound interaction of various transcription factors (TFs). Development and massive application of next generation sequencing technologies for mapping of TF binding sites (BS) in genome, in particular ChIP-seq, provides an opportunity to study gene expression regulation in detail. Widely applied approach of TFBS prediction, the position weight matrix (PWM) relied on a relatively short motifs and proposed the additivity of various positions within potential TFBS [1,2]. Recently, to supplement the verification of potential BSs in ChIP-seq data with PWMs, various approaches that taking into account intra-motifs dependencies were applied for wide-genome application, e.g. BaMM [3] and InMode [4]. These approaches use the markov models (MMs), which neglect the additivity assumption through the concept of the order of markov chain, i.e. a short distance for a given position that may contain other dependent nucleotides. Typically, the order of MM is changed from one to five, that’s way MMs may be referred as to as ‘short-range interaction’ models.
To compare traditional PWMs with BAMM/InMode models we developed the integrated pipeline for ChIP-seq data verification with multiple de novo motif search models.

## Scheme of pipeline

![image](pipeline_scheme.png)

## Requirements

The easiest way to install dependencies is to use [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) or [conda](https://conda.io/projects/conda/en/latest/index.html).
The next command installs all dependencies except BaMM and SiteGA (need to install manually).

1. Create environment and install dependencies
  with _mamba_: `mamba create -n multidena-meme5.4.1 -c bioconda meme=5.4.1 numpy scipy pandas cython r-optparse bioconductor-clusterprofiler bioconductor-org.mm.eg.db bioconductor-org.hs.eg.db bioconductor-org.at.tair.db bioconductor-motifstack bedtools`
  with _conda_: `conda create -n multidena-meme5.4.1 -c bioconda -c conda-forge meme=5.4.1 numpy scipy pandas cython r-optparse bioconductor-clusterprofiler bioconductor-org.mm.eg.db bioconductor-org.hs.eg.db bioconductor-org.at.tair.db bioconductor-motifstack bedtools`

2. Activate environment
 with `mamba`: `mamba activate multidena-meme5.4.1`
 with `conda`: `conda activate multidena-meme5.4.1`

3. Install [BaMM](https://github.com/soedinglab/BaMMmotif2) and [SiteGA](https://github.com/parthian-sterlet/sitega). Also update variable PATH to use BaMM and SiteGA.
4. Clone this repository: `git clone https://github.com/ubercomrade/MultiDeNa.git`
5. Install MultiDeNa
  ```
  cd MultiDeNA/  
  pip3 install -e .
  ```

List of dependencies. if you don't use mamba/conda you can install all dependencies manually

PYTHON:
  * `pip3 install cython numpy scipy pandas`

R:
  * ChIPseek and additional packeges for peaks/scan annotation:
  ```
install.packages("optparse")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.At.tair.db")
  ```
  * motifStack (plot motifs):
  ```
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("motifStack")
  ```

TOOLS:
  * bedtools: https://bedtools.readthedocs.io/en/latest/  version >= 2.26.0

MODELS:
  * BaMM: https://github.com/soedinglab/BaMMmotif2 (BaMMmotif2 is needed to install)
  * ChIPmunk: http://autosome.ru/ChIPMunk/ (ChIPMunk is included in the package)
  * InMoDe: http://jstacs.de/index.php/InMoDe (InMoDe is included in the package)
  * SiteGA: https://github.com/parthian-sterlet/sitega (sitega is needed to install)

OPTIONAL:
  * TomTom: http://meme-suite.org/index.html

## Installation

```  
git clone https://github.com/ubercomrade/MultiDeNa.git  
cd MultiDeNA/  
pip3 install -e .  
```

## Usage
The command `MultiDeNa.py -h` return:

```
MultiDeNa.py -h
usage: MultiDeNa.py [-h] [-b BACKGROUND] [-t TRAIN_SIZE] [-f FPR] [-T TEST_SIZE] [-I INMODE] [-c CHIPMUNK] [-m PATH_TO_MDB]
                    bed N genome output M [M ...]

positional arguments:
  bed                   path to BED file
  N                     Promoters of organism (hg38, mm10, tair10)
  genome                Path to genome fasta file
  output                Output dir
  M                     The comma-separated list of models that will be used in analysis to use.
                        Available options: pwm-chipmunk,pwm-streme,dipwm,bamm,inmode,sitega.
                        At lest two models should be provided. Example: pwm-streme,bamm.

options:
  -h, --help            show this help message and exit
  -b BACKGROUND, --background BACKGROUND
                        Path to background. It is used for de novo and to estimate ROC and PRC. if it is not given background is generated by
                        shuffling
  -t TRAIN_SIZE, --train TRAIN_SIZE
                        Number of peaks for training sample. The default value is 500
  -f FPR, --FPR FPR     FPR, def=1.9*10^(-4)
  -T TEST_SIZE, --test TEST_SIZE
                        Number of peaks for testing sample. It could be any value starting from number of peaks used in traning sample. If
                        parameter is -1 all peaks are used. The default value is -1.
  -l MIN_LENGTH, --min_length MIN_LENGTH
                        Minimal length of motif (default value is 8. Don`t use
                        values less than 6)
  -L MAX_LENGTH, --max_length MAX_LENGTH
                        Maximal length of motif (default value is 20. Don`t
                        use values more than 30)
  -s STEP, --step STEP  The step with which the length of the motif will
                        increase (default value is 4)
  -m PATH_TO_MDB, --motifdatabase PATH_TO_MDB
                        Path to motif database in meme format for TOMTOM. You can get motif database from http://meme-
                        suite.org/doc/download.html
```
### Example run
[Bash script](https://github.com/ubercomrade/MultiDeNA/blob/master/example/multidena_example_run.sh) with example run and tiny data are located in ./example directory. You should run this script in ./example directory.
```
MultiDeNa.py ./foxa2_mm_chr19.bed \
mm10 \
./Mus_musculus.GRCm38.dna.chromosome.19.fa \
./foxa2_mm_chr19 \
pwm-streme bamm \
-t 500 -T -1
```

### Required arguments description

**First positional argument**:
```
bed                   path to BED file
```
You shuld give path to bed file. The file should contain the following columns separated by tab:
1. chromosome
2. start
3. end
4. name
5. score

File can contain additional columns, but they are not used. Information about bed format is avaliable in: https://m.ensembl.org/info/website/upload/bed.html, https://genome.ucsc.edu/FAQ/FAQformat.html#format1

**Second positional argument**:
```
N                     promoters of organism (hg38, mm10, tair10)
```
Value of N can be _hg38_ or _mm10_. It depends on organism used in research

**Third positional argument**:
```
genome                path to genome fasta file
```
Path to genome fasta file.
Reference genomes for mm10/hg38 can be downloaded from https://www.gencodegenes.org/ or https://genome.ucsc.edu/


**Fourth positional argument**:
```
output                output dir
```
Directory path to write results. If directory does not exist it'll be created.


**Fifth positional argument**:
```
N                     The comma-separated list of models that will be used in analysis to use.
                      Available options: pwm-chipmunk,pwm-streme,dipwm,bamm,inmode,sitega.
                      At lest two models should be provided. Example: pwm-streme,bamm.
```
Argument takes different values from the list: pwm-streme, pwm-chipmunk, bamm, inmode, sitega. Several values have to be chosen. Chosen values have to be separated by comma. Example: `pwm-streme,bamm`, `pwm-chipmunk,bamm,inmode`, `pwm-streme,bamm,inmode`, `pwm-streme,sitega`
*IMPORTANT!* Option `bamm` must be used with `pwm-streme` or `pwm-chipmunk`, because BaMM model depends on PWM model. You have to choose more than 1 model in analysis.

### Optional arguments description

```
-t TRAIN_SIZE, --train TRAIN_SIZE
```
The argument `-t/--train` sets the number of peaks that will be used for de novo. The default value is equal to 500.

```
-T TEST_SIZE, --test TEST_SIZE
```
The argument `-T/--test` sets the number of peaks that will be used in analysis. By default all peaks are used.

```
-f FPR, --FPR FPR
```
The argument `-f/--FPR` sets threshold for models to distinguish sites from no-sites. The default value is equal to 0.00019.

```
-l MIN_LENGTH, --min_length MIN_LENGTH
```
Minimal length of motif

```
-L MAX_LENGTH, --max_length MAX_LENGTH
```

Maximal length of motif

```                      
-s STEP, --step STEP
```

The step with which the length of the motif will increase

## Description of results

Multidena will generate output of the following structure (pwm and bamm models are used):
```
├── annotation
│     ├── test
│     │     ├── all_models_GO.tsv
│     │     ├── bamm_annotaion.tsv
│     │     ├── compare_model_GO.tsv
│     │     ├── pwm_annotaion.tsv
│     │     ├── pwm_model_GO.tsv
│     │     └── tools_annotaion.pdf
│     └── train
│           ├── bamm_annotaion.tsv
│           ├── compare_models_GO.tsv
│           ├── pwm_annotaion.tsv
│           ├── pwm_model_GO.tsv
│           └── tools_annotaion.pdf
├── auc
│    ├── bamm
│    │    ├── statistics.txt
│    │    ├── ...
│    └── pwm
│         ├── statistics.txt
│         ├── ...
├── bed
├── fasta
├── models
│     ├── bamm_model
│     │     ├── bamm.hbcp
│     │     ├── bamm_model.ihbcp
│     │     ├── bootstrap_merged.txt
│     │     └── bootstrap.txt
│     ├── pwm_model
│     │     ├── bootstrap_merged.txt
│     │     ├── bootstrap.txt
│     │     ├── pwm_model.meme
│     │     ├── pwm_model.pfm
│     │     ├── pwm_model.pwm
│     │     ├── streme.txt
│     │     └── streme.xml
│     └── thresholds
│         ├── bamm_model_thresholds.txt
│         └── pwm_model_thresholds.txt
├── results
│     ├── combined_scan_test_1.90e-04.pro
│     ├── combined_scan_train_1.90e-04.pro
│     ├── compare_motif_logo.pdf
│     ├── compare_test_1.90e-04_pwm.bamm_counts.tsv
│     ├── compare_train_1.90e-04_pwm.bamm_counts.tsv
│     ├── peaks_classification_test_1.90e-04.tsv
│     └── peaks_classification_train_1.90e-04.tsv
├── scan
│     ├── bamm_test_1.90e-04.bed
│     ├── bamm_train_1.90e-04.bed
│     ├── pwm_test_1.90e-04.bed
│     └── pwm_train_1.90e-04.bed
├── scan-best
│     ├── bamm.scores.txt
│     └── pwm.scores.txt
└── tomtom
      ├── bamm.meme
      ├── bamm.pfm
      ├── bamm.pwm
      ├── bamm.sites.txt
      ├── bamm.tomtom_results
      │    ├── tomtom.tsv
      │    └── tomtom.xml
      ├── pwm.meme
      ├── pwm.pfm
      ├── pwm.pwm
      ├── pwm.sites.txt
      └── pwm.tomtom_results
           ├── tomtom.tsv
           └── tomtom.xml

```
The main results of pipeline are contained in directory `results`.
### `results`

* `combined_scan_{train/test}_{err}` - it's joint profile for all models. It has fasta like format. The file contains results of scan of all models with classification

Example
```
>peaks_0::chr1:164712189-164712712(+) SEQ 1 ## <- first peak
162 170 4.137341427518557 + gataacgg  bamm  1 ## <- start end -log10(err) strand site model_name group (if several models have the same group they are overlapped)
234 242 3.7438654156944446  + gataacag  bamm  2
245 265 3.740915004460945 - cgagcgataaggctgctaac  sitega  3
250 262 3.8460276036525065  - gcgataaggctg  pwm 3
252 260 3.9555972318346377  - gataaggc  bamm  3
357 377 4.715008011604894 - cttccgatatcaaccgggaa  sitega  4
440 460 4.364457279299288 - atggtgattaacatccccct  sitega  5
>peaks_1::chr3:93470273-93470886(+) SEQ 2 ## <- second peak
106 126 3.9101503930828088  + actttgatatttcatgtaga  sitega  1
>peaks_2::chr17:36862394-36862756(+)  SEQ 3
129 137 3.7257217215553267  + gataacac  bamm  1
193 205 3.8499146293487376  + cagataaccgag  pwm 2
195 203 3.775843524673776 + gataaccg  bamm  2
```

* `compare_motif_logo.pdf` - the file contains logos for all models

* `compare_{train/test}_{err}_{model1}.{model2}_counts.tsv` - the file contains pair comparison of models. Classification with taking into account overlapping
Example:
```
pwm: against bamm bamm: against pwm overlapped:pwm,bamm not_overlapped:pwm,bamm no_sites:pwm,bamm number_of_peaks
328               301               2123                48                      4360              7160
```

* `peaks_classification_*`- the file contains classification of all peaks for all model. Classification without taking into account overlapping
Example (three model):
```
pwm bamm  sitega  bamm_and_pwm  pwm_and_sitega  bamm_and_sitega bamm_and_pwm_and_sitega number_of_peaks
1793  2494  4956  6403  459 709 1716  46458
```

## Useful links

 * [Cistrome](http://cistrome.org/ap/) [5]

## Reference

[1]	Benos P.V. et al. (2002) Additivity in protein-DNA interactions: how good an approximation is it? Nucleic Acids Res., 30(20):4442-4451.  
[2]	Srivastava D and Mahony S. (in press) Sequence and chromatin determinants of transcription factor binding and the establishment of cell type-specific binding patterns. Biochim Biophys Acta Gene Regul Mech., 194443. doi: 10.1016/j.bbagrm.2019.194443.  
[3]	Siebert M. and Söding J. (2016) Bayesian Markov models consistently outperform PWMs at predicting motifs in nucleotide sequences. Nucleic Acids Res., 44(13):6055–6069.  
[4]	Eggeling R. et al. (2017) InMoDe: tools for learning and visualizing intra-motif dependencies of DNA binding sites. Bioinformatics, 33(4):580-582.  
[5]	Mei S. et al. (2017) Cistrome Data Browser: a data portal for ChIP-Seq and chromatin accessibility data in human and mouse, Nucleic Acids Res., 45(D1):, D658–D662,.  

## License

Copyright (c) 2020 Anton Tsukanov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
