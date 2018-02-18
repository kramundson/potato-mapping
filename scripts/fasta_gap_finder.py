#!/usr/bin/env python3

# USAGE: python fasta_gap_finder.py -f fastafile -n buffer

# Import necessary packages
import argparse
import re
from Bio import SeqIO

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-f", dest="fasta")
parser.add_argument("-n", dest="buffer", type=int, default=0)
args = parser.parse_args()

# Open FASTA, search for masked regions ± buffer, print in BED3 format
# If wanting to keep non-gap region, use the original FASTA and BED file of this script in bedtools subtract
with open(args.fasta) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        for match in re.finditer('N+', str(record.seq)):
            print(record.id, match.start()-args.buffer, match.end()+args.buffer, sep='\t')
