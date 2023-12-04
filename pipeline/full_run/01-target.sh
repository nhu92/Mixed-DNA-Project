#!/bin/bash
#SBATCH -J target_run
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
threads=64
read1=.R1.fastq.gz
read2=.R2.fastq.gz
mega353=mega353.fasta
proj_name=
ref_alignment=ref # The directory of reference alignment

# ---

fastp -i ${read1} -I ${read2} -o ${read1}.trimmed.fastq.gz -O ${read2}.trimmed.fastq.gz -j fastp.json -h fastp.html
mkdir hyb_output
hybpiper assemble -t_dna ${mega353} \
        -r ${read1}.trimmed.fastq.gz ${read2}.trimmed.fastq.gz \
        --prefix ${proj_name} \
        --bwa \
        --cpu ${threads} \
        -o ./hyb_output

find hyb_output/${proj_name} -maxdepth 1 -type d -exec basename {} \; | grep -Po "\d+" > gene_list.txt
mkdir ./exon_extracted
while read line; do python split_exon_extract.py ./hyb_output/${proj_name} $line ./exon_extracted 0.8; done < gene_list.txt

input_exon=./exon_extracted

mkdir ./phylo_results
# A batch run to generate trees for all genes:
while read gene_name_shorter
do
        ls ${input_exon}/*${gene_name_shorter}*.fasta > ./exon.list.txt
        i=1
        while read exons
        do
                # MAFFT alignment and trim over-gapped
                mafft --preservecase --maxiterate 1000 --localpair --adjustdirection --thread ${threads} --addfragments ${exons} ${ref_alignment}/${gene_name_shorter}.fasta > phylo_results/${gene_name_shorter}_exon_${i}_aligned.fasta
                trimal -in phylo_results/${gene_name_shorter}_exon_${i}_aligned.fasta -out phylo_results/${gene_name_shorter}_exon_${i}_trimmed.fasta -gt 0.5

                # Tree construction and comparison
                fasttree -gtr -gamma -nt phylo_results/${gene_name_shorter}_exon_${i}_trimmed.fasta > phylo_results/${gene_name_shorter}_exon_${i}.tre
                ((i=i+1))
        done < exon.list.txt
done < gene_list.txt

rm exon.list.txt
# Distance Calc (in progress)
mkdir all_trees

cp phylo_results/*.tre ./all_trees

tree_dir=all_trees
while read gene_name_shorter
do
        ls "${tree_dir}/${gene_name_shorter}"*"tre" > ./loop.treelist.txt
        i=1
        while read filename
        do
                python matrix.py -t ${filename} -o ./all_trees/${gene_name_shorter}.${i}.matrix
		cp ./all_trees/${gene_name_shorter}.${i}.matrix ./all_trees/${gene_name_shorter}.${i}.cleaned.csv
		((i=i+1))
        done < loop.treelist.txt
done < ./gene_list.txt

rm loop.treelist.txt

python dist2Z.py all_trees ./${proj_name}.summary_dist.csv ${proj_name}
python group_sum.py ./${proj_name}.summary_dist.csv ./${proj_name}.cumulative_dist.csv
grep -v ${proj_name} ${proj_name}.cumulative_dist.csv | cut -d ',' -f2 > ${proj_name}.Zaxis.csv

paste --delimiter=',' ${proj_name}_refPCoA.csv ${proj_name}.Zaxis.csv > ${proj_name}.3d.csv
python contour_optima_normalized.py ${proj_name}.3d.csv ${proj_name}.contour.norm.svg 0.5
paste --delimiter=',' ${proj_name}_refPCoA23.csv ${proj_name}.Zaxis.csv > ${proj_name}.3d23.csv
python contour_optima_normalized.py ${proj_name}.3d23.csv ${proj_name}.contour23.norm.svg 0.5
