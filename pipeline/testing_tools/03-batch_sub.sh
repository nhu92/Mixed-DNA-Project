#!/bin/bash


# Change the 1 and 16 till we finish 100 submissions.
i=1
while read line
do
    if [ $i -ge 1 ] && [ $i -le 16 ]
    then
        mkdir -p /lustre/scratch/nhu/202312/3ndpara_Test/test_lab/"${line}"/log
        cd /lustre/scratch/nhu/202312/3ndpara_Test/test_lab/"${line}"
	sbatch 01-target.sh
    fi
    i=$((i+1))
done < mix_namelist.txt

