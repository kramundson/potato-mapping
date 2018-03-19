import sys

"""
USAGE:

To be run in snakemake pipeline from ./potato-mapping in rule

script:
    scripts/make_intervals.py <chunksize>
"""

def make_intervals(chunksize):

    """
    Reading from file potato_dm_v404_all_pm_un.dict, return 1 line for each genomic
    window of fixed size. The size is determined from the argument chunksize. 
    Each output line is a 1-based interval of a specific chromosome, e.g.,
    
    chr01:1-1000000 # first interval line
    chr01:1000001-2000000 # second interval line
    
    For the final interval of a chromosome, return a shortened interval.
    If the length of chr01 is 88663952:
    
    chr01:88000000-88663952 # final interval file of chromosome 1
    """

    chunk_int = int(chunksize)
    ofh="potato_"+sys.argv[1]+".intervals"
    o=open(ofh, 'w')
    with open("data/genome/potato_dm_v404_all_pm_un.dict") as f:
        for line in f:
            if line.startswith("@SQ\tSN:"):
                l = line.split('\t') # consider reducing this to one line
                chrom = l[1].split(":")[1]
                end = int(l[2].split(":")[1])
                for i in range(1, end, chunk_int):
                    if i + chunk_int > end:
                        outstr = chrom + ":" + str(i) + "-" + str(end) + "\n"
                    else:
                        outstr = chrom + ":" + str(i) + "-" + str(i+chunk_int-1) + "\n"
                    o.write(outstr)
    o.close()

if __name__ == "__main__":
    make_intervals(sys.argv[1])
