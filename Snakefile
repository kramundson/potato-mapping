import re
import pandas as pd
shell.executable("bash")

configfile: "config.yaml"

units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels]) # enforce str in index

def is_single_end(sample,unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def get_fastq(wildcards):
    return "data/reads/"+units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample should be aligned as such
        return expand("data/trimmed/{sample}-{unit}-{group}.fastq.gz",
            group=[1,2], **wildcards)
    # single end sample
    return "data/trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)

rule all:
    input:
        config["genome"],
        config["genome"]+".bwt",
        "data/intervals/gatk-haplocaller.intervals",
        config["vcfs"]["diploid"],
        config["vcfs"]["tetraploid"],
        config["vcfs"]["combined"]

# TODO: make sure all thread usage in all rules is specified in config.yaml

include: "rules/get_genome.rules"
include: "rules/get_SRA_reads.rules"
include: "rules/cutadapt.rules"
include: "rules/cutadapt_pe.rules"
include: "rules/bwa_index.rules"
include: "rules/align.rules"
include: "rules/mark_duplicates.rules"
include: "rules/samtools_index.rules"
include: "rules/samtools_faidx.rules"
include: "rules/gatk4_fasta_dict.rules"
include: "rules/make_intervals.rules"
include: "rules/gatk4_parallel_fofns.rules"
include: "rules/gatk4_haplotypecaller_diploid_parallel.rules"
include: "rules/gatk4_haplotypecaller_tetraploid_parallel.rules"
include: "rules/gatk4_combine_sample_gvcfs.rules"
include: "rules/gatk4_genotype_gvcfs.rules"
