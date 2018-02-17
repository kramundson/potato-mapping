#!/bin/bash

cut -f 2 units.tsv | grep "^SRR" | xargs fastq-dump --gzip -B --split-3 -O ./data/reads

