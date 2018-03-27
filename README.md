Potato mapping and variant calling from short read sequencing

Before running, install the Python3 version of Miniconda

[Installation instructions](http://snakemake.readthedocs.io/en/stable/tutorial/setup.html)

Once Miniconda is installed, install software dependencies using included environment.yaml file.
To set up a conda environment, see the tutorial 

[Environment creating instructions](https://conda.io/docs/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)

To activate this environment on the command line:

```
source activate <name_of_environment>
```

This Snakemake workflow runs the following:

1. Download of DM1-3 reference genome from http://solanaceae.plantbiology.msu.edu/pgsc_download.shtml
2. Download Illumina sequence data for all publicly available samples in units.tsv (those that start with "SRR")
3. Read quality and adapter trimming
4. BWA mem alignment to reference, including sam to bam conversion, bam sort, and samtools index.
    Paired end reads should be uninterleaved and are are parsed from units.tsv by having
    the appropriate file name in column fq2 of units.tsv. See units.tsv for example.
    For single end reads, set the fq2 field of that row to "NaN"
5. HaplotypeCaller from GATK4 in GVCF mode on each sample, scattered across intervals.
6. CombineGVCFs from biological samples, using GATK4 CombineGVCFs
7. Joint genotype calling on per-sample GVCFs using GATK4 GenotypeGVCFs

Outputs:
1. One GVCF file per sample and its index
2. Population VCF file
3. Currently, all intermediate files (e.g., trimmed reads, unprocessed bams, region GVCFs) 
    are also kept. Uses too much disk space, so this will likely go away soon.

Configuration:
1. Modify units.tsv to suit your needs. Each column specifies:
    sample: unique biological sample
    unit: unique combination of biological sample, library prep, and sequencing run
    fq1: name of forward read
    fq2: name of reverse read (enter NaN here if reads are single-ended)
    parhap: not actually used, yet
    
    Note: Avoid using the "-" character in the sample and unit fields.

2. Place reads fq1 and fq2 in the subfolder ```data/reads/```
3. Modify parameters, thread usage, and names of target output files in config.yaml
4. Snakemake will automatically spawn jobs when running on a cluster. If desired, you can
change the memory and CPU requirements of each job (as well as other params) by
modifying the file cluster.yaml
5. Run pipeline. In a cluster, the job can be submitted with the following command:

```sbatch runSnakes.slurm```

This command will run two Snakefiles in succession. The first, init_genome.snakes, 
downloads the reference, generates reference index files, and sets up intervals that GATK4
will operate on in parallel. The second file, Snakefile, downloads and processes sequencing
reads.