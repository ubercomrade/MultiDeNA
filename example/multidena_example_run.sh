#!/usr/bin/env bash

if [[ -e Mus_musculus.GRCm38.dna.chromosome.19.fa.gz ]]
then
    gunzip Mus_musculus.GRCm38.dna.chromosome.19.fa.gz
fi

MultiDeNa.py ./foxa2_mm_chr19.bed \
mm10 \
./Mus_musculus.GRCm38.dna.chromosome.19.fa \
./foxa2_mm_chr19 \
pwm-streme bamm \
-t 500 -T -1
