#!/bin/bash
#
#SBATCH --job-name=snake_mapping
#SBATCH -c 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=4G # Memory pool for all cores in MB (see also --mem-per-cpu)
#SBATCH --time 3-00:00:00 # 3 day runtime
#SBATCH -p gc # Partition to submit to
#SBATCH -o 2018_0321_snake.out # File to which STDOUT will be written
#SBATCH -e 2018_0321_snake.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=kramundson@ucdavis.edu # Email to which notifications will be sent

date
hostname
cd /share/comailab/kramundson/potato-mapping # modify path as needed

source activate gatk-potato # alter as needed

snakemake -j 7 -s init_genome.snakes --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time} -c {cluster.c} --mem-per-cpu {cluster.mempercpu}" -k -w 120
snakemake -j 999 --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time} -c {cluster.c} --mem-per-cpu {cluster.mempercpu}" -k -w 120

source deactivate gatk-potato

date
