import sys
from Bio import SeqIO

"""
USAGE: script.py fastafile > cleanfastafile
"""

with open(sys.argv[1], 'r') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        print(">{}\n{}".format(record.id,record.seq))