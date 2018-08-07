# Run init_genome.snakes before running this Snakefile
# This Snakefile handles dataset-specific analysis, assuming that reference genome
# file dependencies and file of output filenames have been made using init_genome.snakes

import os, re
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
        # paired-end sample
        return expand("data/trimmed/{sample}-{unit}-{group}.fastq.gz",
            group=[1,2], **wildcards)
    # single end sample
    return "data/trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)

def get_vcfs(wildcards):
    intervals = [line.rstrip('\n') for line in open("data/intervals/gatk-haplocaller.intervals", 'r')]
    vcfs = ["-I=data/calls/gatk/vcf_by_region/all-samples-{}.vcf".format(x) for x in intervals]
    return vcfs

rule all:
    input:
        "data/calls/all-calls.vcf"

include: "rules/gatk4_merge_vcfs.rules"
include: "rules/gatk4_genotype_region_gvcfs.rules"
include: "rules/gatk4_merge_gvcfs_by_region.rules"
include: "rules/gatk4_haplotypecaller_tetraploid_cluster.rules"
include: "rules/gatk4_haplotypecaller_diploid_cluster.rules"
include: "rules/samtools_index.rules"
include: "rules/samtools_merge.rules"
include: "rules/mark_duplicates.rules"
include: "rules/merge_aligned_pear.rules"
include: "rules/align_pear_se.rules"
include: "rules/align_pear_pe.rules"
include: "rules/pear.rules"
include: "rules/cutadapt_pe.rules"
include: "rules/cutadapt.rules"
include: "rules/get_SRA_reads.rules"
