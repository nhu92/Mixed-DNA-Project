#!/bin/bash
#SBATCH -J 100prep
#SBATCH -p nocona
#SBATCH -o log/%x.out
#SBATCH -e log/%x.err
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=3G

# Constants. Please edit them to match the data structure of your own directory. Don't include the last "/" in the path.

# A batch run for all genes:
while read mix
do
	mkdir ../test_lab/${mix}
	cp -r ../full_run/* ../test_lab/${mix}
	sed -i s/target_run/${mix}/g ../test_lab/${mix}/01-target.sh
	sed -i s/.R1.fastq.gz/${mix}.R1.fastq.gz/g ../test_lab/${mix}/01-target.sh
	sed -i s/.R2.fastq.gz/${mix}.R2.fastq.gz/g ../test_lab/${mix}/01-target.sh
	sed -i s/proj_name=/proj_name=${mix}/g ../test_lab/${mix}/01-target.sh
done < mix_namelist.txt


