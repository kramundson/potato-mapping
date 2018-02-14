#!/bin/bash

# subsample to get ~0.5x coverage with 2x150nt paired-end reads
# for the potato genome, this amounts to 1.4 million pairs
module load seqtk/1.0-r63-dirty

cd ~/LOP_highcov/data/trimmed/
for i in LOP868-SRR6123032-{1,2}.fastq.gz
do
    BASE=$(basename $i '.gz')
    echo $BASE
    seqtk sample -s100 $i 1400000 > 'sub_seed100_'$BASE
done
