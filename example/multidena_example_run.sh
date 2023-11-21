#!/usr/bin/env bash

MultiDeNa.py ./PEAKS042778_PIF5_chr_1.bed \
tair10 \
./Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa \
./PEAKS042778_PIF5_results \
pwm-streme bamm \
-t 500 -T -1
