#!/usr/bin/python3
# tarean_headers.py
# 31 January 2018
# Kirk Amundson

"""
Reformat NCBI SRA fastq header for TAREAN analysis.
File should be interleaved paired-end sequencing reads in FASTA format
"""

# USAGE
# script.py input.fasta

import sys

fh = "tarean_"+sys.argv[1]
print(fh)
o = open(fh, 'w')

# Open and read 4 lines at a time
with open(sys.argv[1]) as f:
    while 1:
        head1 = f.readline().split()
        if head1 == []:
            break
        else:
            seq1= f.readline()
            head2 = f.readline().split()
            seq2 = f.readline()
            head1out = head1[0] + '/1' + '\n'
            head2out = head2[0] + '/2' + '\n'
            o.write(head1out+seq1+head2out+seq2)

o.close()