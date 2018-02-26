import pandas as pd

def get_diploids():
    diploids = units.loc[units["parhap"] == "haploid"].index.values
    return ["data/dedup/" + d[0] + "-" + d[1] + "-dedup.sorted.aln.bam" for d in diploids]
    
def get_tetraploids():
    tetraploids = units.loc[units["parhap"] == "mother"].index.values
    return ["data/dedup/" + t[0] + "-" + t[1] + "-dedup.sorted.aln.bam" for t in tetraploids]
    
# def get_diploid_bams(wildcards):
#     return "data/reads/" + units.loc[(wildcards.sample, wildcards.unit), ["ploidy"]==2] + "-dedup-sorted.aln.bam"
# 
# def get_tetraploid_bams(wildcards):
#     return "data/reads/" + units.loc[(wildcards.sample, wildcards.unit), ["ploidy"]==4] + "-dedup-sorted.aln.bam"

# def get_diploid_bams():
#     return "data/reads/" + units.loc[units["parhap"]=="haploid"] + "-dedup-sorted.aln.bam"
#     
# def get_tetraploid_bams():
#     return "data/reads/" + units.loc[units["parhap"]=="mother"] + "-dedup-sorted.aln.bam"

def is_single_end(sample,unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def is_diploid(sample,unit):
    return units.loc[(sample,unit), units["ploidy"] == 2]
    
def is_tetraploid(sample,unit):
    return units.loc[(sample,unit), units["ploidy"] == 4]

def get_tetraploids(wildcards):
    if is_tetraploid(**wildcards):
        return "data/dedup/{sample}-{unit}-dedup-sorted.aln.bam".format(**wildcards)

def get_diploids(wildcards):
    if is_diploid(**wildcards):
        return "data/dedup/{sample}-{unit}-dedup-sorted.aln.bam".format(**wildcards)


units = pd.read_table("units_tester.tsv", index_col=["sample","unit"], dtype=str)
units


print(is_single_end)

# print(get_diploids())
# print(get_tetraploids())
# print(get_diploid_bams())
# print(get_tetraploid_bams())
