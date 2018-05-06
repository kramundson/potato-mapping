# Run init_genome.snakes before running this Snakefile
# This Snakefile handles dataset-specific analysis, assuming that reference genome
# file dependencies and file of output filenames have been made using init_genome.snakes

import re
import pandas as pd
shell.executable("bash")

configfile: "config.yaml"

units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

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

def get_region_gvcfs(wildcards):
    return [line.rstrip('\n') for line in open("fofn/{sample}-{unit}.fofn".format(**wildcards))]

rule all:
    input:
        config["vcfs"]["diploid"],
        config["vcfs"]["tetraploid"],
        config["vcfs"]["combined"]

include: "rules/gatk4_genotype_gvcfs.rules"
include: "rules/gatk4_combine_sample_gvcfs.rules"
include: "rules/gatk4_gather_region_gvcfs.rules"
include: "rules/gatk4_haplotypecaller_tetraploid_cluster.rules"
include: "rules/gatk4_haplotypecaller_diploid_cluster.rules"
include: "rules/samtools_index.rules"
include: "rules/mark_duplicates.rules"
include: "rules/align.rules"
include: "rules/cutadapt_pe.rules"
include: "rules/cutadapt.rules"
include: "rules/get_SRA_reads.rules"
