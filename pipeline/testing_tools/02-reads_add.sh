#!/bin/bash
#SBATCH -J 100reads
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
	echo "Value of mix: ${mix}"
	cp ../reads_mix/${mix}* ../test_lab/"${mix}"/
done < mix_namelist.txt


