#!/bin/bash
#SBATCH -J _rescue
#SBATCH -p nocona
#SBATCH -o log/%x.out
#SBATCH -e log/%x.err
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=3G

# This is a pipeline for analyzing raw reads from a mixed sample to the identification.
# This pipeline is to treat with the mixed reads. There is another pipeline to generate the reference panel.

# Declare the constants
proj_name=

# ---

python dist2Z.py all_trees ./${proj_name}.summary_dist.csv ${proj_name}
python group_sum.py ./${proj_name}.summary_dist.csv ./${proj_name}.cumulative_dist.csv
grep -v ${proj_name} ${proj_name}.cumulative_dist.csv | cut -d ',' -f2 > ${proj_name}.Zaxis.csv

