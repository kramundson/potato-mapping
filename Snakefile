import pandas as pd
shell.executable("bash")

configfile: "config.yaml"

samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels]) # enforce str in index

def is_single_end(sample,unit): # used to specify single vs. paired end during alignment
    return pd.isnull(units.loc[(sample, unit), "fq2"]) # fq2 comes from units.tsv

def get_fastq(wildcards):
    return "data/reads/"+units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def get_trimmed(wildcards): # keep here, move to align.smk once made
    if not is_single_end(**wildcards):
        # paired-end sample should be aligned as such
        return expand("data/trimmed/{sample}-{unit}-{group}.fastq.gz",
            group=[1,2], **wildcards)
    # single end sample
    return "data/trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)

#######an earnest attempt at integrating fastq dump of SRA data into this pipeline#######
# samples = pd.read_table(config["samples"], index_col="sample")
# units = pd.read_table(config["units"], index_col="sample", dtype=str)
# units.index = units.index.set_levels([i.astype(str) for i in units.index.levels]) # enforce str in index

# def is_single_end(sample):
#     return pd.isnull(units.loc[(sample), "fq2"]) # fq2 comes from units.tsv
# 
# def get_fastq(wildcards):
#     return "data/reads/"+units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
# 
# def get_trimmed(wildcards):
#     if not is_single_end(**wildcards):
#         # paired end sample
#         return expand("data/trimmed/{sample}-{unit}--{group}.fastq.gz",
#             group=[1,2], **wildcards)
#     # single end sample
#     return "data/trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)
######End earnest attempt, resume normal pipeline######

rule all: # TODO clean this up
    input:
        config["genome"],
        "data/genome/potato_dm_v404_all_pm_un.fasta.bwt",
        "data/dedup/LOP868_004-KFRAG_00003H-dedup-sorted.aln.bam.bai",
        "data/dedup/LOP868_064-KFRAG_00028H-dedup-sorted.aln.bam.bai",
        "data/dedup/LOP868_305-KFRAG_00092H-dedup-sorted.aln.bam.bai",
        "data/dedup/LOP868-SRR6123032-dedup-sorted.aln.bam.bai",
        "data/dedup/IVP101-SRR6123183-dedup-sorted.aln.bam.bai",
        "data/dedup/PL4-SRR6123031-dedup-sorted.aln.bam.bai",
        "data/calls/all-calls.vcf"

rule get_potato_genome:
    output:
        "data/genome/potato_dm_v404_all_pm_un.fasta"
    shell:
        "./scripts/get_potato_genome.sh"

#rule datagrab: # todo get this running else omit
#    input:
#        get_SRRid
#    output:
#        "data/reads/SRR6123031_1.fastq.gz",
#        "data/reads/SRR6123031_2.fastq.gz",
#        "data/reads/SRR6123032_1.fastq.gz",
#        "data/reads/SRR6123032_2.fastq.gz",
#        "data/reads/SRR6123183_1.fastq.gz",
#        "data/reads/SRR6123183_2.fastq.gz"
#    shell:
#        "fastq-dump -B -I --gzip --split-3 -O ./data/reads {input}"

rule cutadapt_pe:
    input:
        get_fastq
    output:
        fastq1="data/trimmed/{sample}-{unit}-1.fastq.gz",
        fastq2="data/trimmed/{sample}-{unit}-2.fastq.gz",
        qc="data/trimmed/{sample}-{unit}.qc.txt"
    threads: 12
    params:
        "-a {} -A {} {}".format(config["adapter"], config["adapter"], config["params"]["cutadapt-pe"])
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.17.4/bio/cutadapt/pe"

rule cutadapt:
    input:
        get_fastq
    output:
        fastq="data/trimmed/{sample}-{unit}.fastq.gz",
        qc="data/trimmed/{sample}-{unit}.qc.txt"
    threads: 12
    params:
        "-a {} {}".format(config["adapter"], config["params"]["cutadapt-se"])
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.17.4/bio/cutadapt/se"

rule bwa_index: # need to add this rule and run properly
    input:
        "data/genome/potato_dm_v404_all_pm_un.fasta"
    output:
        "data/genome/potato_dm_v404_all_pm_un.fasta.bwt"
    shell:
        "bwa index {input}"

rule align:
    input:
        reads=get_trimmed,
        ref=config["genome"],
        index="data/genome/potato_dm_v404_all_pm_un.fasta.bwt"
    output:
        "data/aligned_reads/{sample}-{unit}-sorted.aln.bam"
    log:
        "logs/bwa_mem/{sample}-{unit}.log"
    params:
        rg="'@RG\\tID:{unit}\\tSM:{unit}'"
    shell: # todo optimize CPU and memory usage
        "bwa mem -R {params.rg} -t 12 {input.ref} {input.reads} | "
        "samtools sort -@4 -m 4G -o {output} -"

# uses snakemake wrapper repository, https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/picard/markduplicates.html
# Runs out of java heap space if using miniconda default.
# To correct, add desired memory to jvm_mem_opts in the file:
# miniconda3/envs/potato/share/picard-2.14.1-0 
# I doubled the heap space from the default: -Xms512m -Xmx1g. This seemed to work fine.

rule mark_duplicates:
    input:
        "data/aligned_reads/{sample}-{unit}-sorted.aln.bam"
    output:
        bam="data/dedup/{sample}-{unit}-dedup-sorted.aln.bam",
        metrics="data/dedup/{sample}-{unit}-metrics.txt"
    log:
        "logs/picard/{sample}-{unit}.log"
    params:
        "REMOVE_DUPLICATES=true",
        "TMP_DIR=./temp",
        "ASSUME_SORT_ORDER=coordinate",
        "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000"
    wrapper:
        "0.19.3/bio/picard/markduplicates"

rule samtools_index:
    input:
        "data/dedup/{sample}-{unit}-dedup-sorted.aln.bam"
    output:
        "data/dedup/{sample}-{unit}-dedup-sorted.aln.bam.bai"
    shell:
        "samtools index {input}"

rule freebayes: # todo implement snakemake wrapper
    input:
        ref=config["genome"]
    output:
        "data/calls/all-calls.vcf"
    log:
        "logs/freebayes/all-calls.log"
    shell:
         "freebayes -f {input.ref} data/dedup/*.bam > {output}"
         
rule bedtools_init:
    input:
        ref=config["genome"]
    output:
        "data/genome/potato_dm_v404_all_pm_un_1Mb_windows.bed",
        "data/genome/potato_dm_v404_all_pm_un.genome"
    shell:
        "./scripts/bedtools-init.sh"

rule bedtools_coverage:
    input:
        bam="data/dedup/{sample}-{unit}-dedup-sorted.aln.bam",
        genome=config["genome"],
        windows="data/genome/potato_dm_v404_all_pm_un_1Mb_windows.bed"
    output:
        "data/dedup/{sample}-{unit}_coverage.bed"
    shell:
        "bedtools coverage -sorted -nonamecheck -header -g {input.genome} -a {input.windows} -b {input.bam} > {output}"

rule bedfiles:
    input:
    output:
        "data/dedup/bedfiles.txt"
    shell:
        "ls data/dedup/*_coverage.bed > data/dedup/bedfiles.txt"

rule dosage_plot: # todo: a more snakemake way to do this is make the script consider case and control only, rather than make all at once
    input:
        control=config["controlfile"],
        bedfiles="data/dedup/bedfiles.txt"
    output:
        "data/plots/{sample}-{unit}_dosage_plot.pdf"
    shell:
        "Rscript scripts/dosage_plots.R {input.bedfiles} {input.control}"