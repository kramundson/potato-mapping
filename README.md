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

