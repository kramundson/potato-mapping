#!/bin/bash
# fastInterleave.sh
# 31 January 2018
# Interleave a pair of uninterleaved fastq.gz files

# USAGE:
# fastInterleave.sh reads-1.fastq.gz reads-2.fastq-gz

echo $1
echo $2

BASE=$(basename $1 '-1.fastq.gz')
echo $BASE'.fastq.gz'

paste <(gunzip -c $1 | paste - - - -) <(gunzip -c $2 | paste - - - -) | tr '\t' '\n' | gzip -c > $BASE'.fastq.gz'
