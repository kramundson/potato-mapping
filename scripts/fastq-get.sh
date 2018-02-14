#!/bin/bash
module load sratoolkit/2.8.1-3

for i in {SRR6123183,SRR6123032,SRR6123031}
do
    fastq-dump --gzip -B --split-3 -O ../data/reads/ $i
done
