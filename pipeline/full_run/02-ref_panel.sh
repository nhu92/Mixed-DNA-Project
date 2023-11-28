#!/bin/bash
#SBATCH -J ref_panel
#SBATCH -p nocona
#SBATCH -o log/%x.out
#SBATCH -e log/%x.err
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=3G

# This is a pipeline for transforming reference alignment to a PCoA panel.
# This pipeline is to treat with the reference panel. There is another pipeline to generate the mixed input distance matrix.


# Constants. Please edit them to match the data structure of your own directory. Don't include the last "/" in the path.
iqtree_dir=~/software/iqtree/iqtree-2.2.2.7-Linux/bin
ref_alignment=ref # The directory of reference alignment
proj_name=

# A batch run for all genes:
ls ${ref_alignment} > ./allgenes.txt
mkdir ref_tree
while read gene
do
        # Tree construction and comparison
        ${iqtree_dir}/iqtree2 -s ${ref_alignment}/${gene} -m MFP -bb 1000 -redo -pre ref_tree/${gene}
        nw_ed ref_tree/${gene}.treefile 'i & b<50' o > ref_tree/${gene}.collapsed.tre
done < ./allgenes.txt

cat ref_tree/*.collapsed.tre | sed 's/_R_//g' | sed 's/_\([0-9]\{4\}\)//g' > ${proj_name}_merged.collapsed.tre
python mean_distance_matrix.py ${proj_name}_merged.collapsed.tre ${proj_name}_meandist.csv
python dist2PCoA.py ${proj_name}_meandist.csv ${proj_name}_refPCoA.csv ${proj_name}_refPCoA.svg
python dist2PCoA23.py ${proj_name}_meandist.csv ${proj_name}_refPCoA23.csv ${proj_name}_refPCoA23.svg