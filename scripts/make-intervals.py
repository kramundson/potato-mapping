def scatter_gvcf(): # remove wildcards arg if debugging without Snakemake
# def scatter_gvcf(wildcards): # note, will need to add wildcards arg when using with Snakemake
    chunk_int=int(1.5e7) # temp hardcode at 15Mb
    intervals=[] # temp initiate empty list. Populate and pass to return.
    with open("../data/genome/potato_dm_v404_all_pm_un.dict") as f:
        for line in f:
            if line.startswith("@SQ\tSN:"):
                l = line.split('\t')
                chrom=l[1].split(":")[1]
                end=int(l[2].split(":")[1])
                for i in range(0, end, chunk_int):
                    if i + chunk_int > end:
                        out= chrom + ":" + str(i) + "-" + str(end)
                    else:
                        out=chrom + ":" + str(i) + "-" + str(i+chunk_int)
                    intervals.append(out)
#     return intervals
    return expand("path/gatk/2x_{sample}-{unit}-{window}.g.vcf", window=[intervals], **wildcards)
    # make a list of all genomic regions. I've done this in the past
    #return expand("data/calls/gatk/2x_{sample}-{unit}-{region}.g.vcf", region=

if __name__ == "__main__":
    print(scatter_gvcf())