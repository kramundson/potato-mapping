Snakemake pipeline for read alignment to potato reference genome, FreeBayes variant calling, and
dosage plot generation. Dosage plots intended for use with low-coverage sequencing data.

If miniconda3 is not installed, follow the instructions at this URL:
http://snakemake.readthedocs.io/en/stable/tutorial/setup.html

Once miniconda is installed, you can  build a suitable environment from environment .yaml

To run this analysis, place your reads in data/reads/ and update both samples.tsv and units.tsv
You should provide unique unit identifiers for your reads. Multiple units may come from the same sample.

In other words, consider unit a unique combination of biological sample, sequencing library, and sequencing run.
Consider a sample a unique biologcical sample that may be represented by multiple libraries and sequencing runs.
Empty files with the appropriate format are located in data/reads that should enable a Snakemake dry run to be performed.

Notes for running GATK 3.x:
Due to licensing restrictions, the bioconda recipe of GATK 3.8 does not work right out of the box.

To install, first set up the conda environment, then download GATK 3.8.
Had to download over web browser, which yielded the following .bz2:

Name: GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
md5sum: 8cfe4517dad2bea5b79c46ed468d62dc
Upload to your working environment with scp or equivalent

To install:

```
tar -xjvf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
gatk-register /path/to/GenomeAnalysisTK-3.8-1-0-gfc1c3ef/GenomeAnalysisTK.jar
```

I uploaded the .bz2 to ./scripts, making my version of the command

```
scp path/on/local/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 path/on/server/potato-mapping/scripts
tar -xjvf ./scripts/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
gatk-register GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
```