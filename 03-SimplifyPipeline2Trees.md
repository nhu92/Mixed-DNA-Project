# An easier way to filter the data. Use HybPiper again!

Utilize the exonerate data to extract conservative region across all the hits contigs.
---

As the original pipeline works, the remaining issue now is to expand the "in" pool of gene tree construction. Rather than adjust the filtering criteria, the more valid way is to direct extract the targeted gene regions from the output tsv from the exonerate run of HybPiper.

Let's check the format of the exonerate output. The file is stored in `{hybpiper_dir}/{gene_name}/{hybpiper_proj_name}/exonerate_stats.tsv`.

Here is a preview of this file:

|                                       | query_id  | query_length | hit_id                            | query_HSP_range_limits_original | query_HSP_range_limits_trimmed | query_HSPFragment_ranges                                                       | hit_percent_similarity_original | hit_percent_similarity_trimmed | hit_strand | hit_HSP_range_limits_original | hit_HSP_range_limits_trimmed | hit_HSPFragment_ranges_original                                                              | hit_HSPFragment_ranges_trimmed                                                               | 3-prime_bases_trimmed |
|---------------------------------------|-----------|--------------|-----------------------------------|---------------------------------|--------------------------------|--------------------------------------------------------------------------------|---------------------------------|--------------------------------|------------|-------------------------------|------------------------------|----------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|-----------------------|
| Hits filtered > 55 percent similarity | VLNB-5974 | 234          | NODE_8_length_1347_cov_20.507377  | (0, 132)                        | (0, 132)                       | [(0, 5), (5, 21), (21, 34), (35, 68), (68, 132)]                               | 92.48                           | 92.48                          | 1          | (36, 1303)                    | (36, 1303)                   | [(36, 51), (290, 338), (468, 508), (603, 704), (1111, 1303)]                                 | [(36, 51), (290, 338), (468, 508), (603, 704), (1111, 1303)]                                 | N/A                   |
|                                       | VLNB-5974 | 234          | NODE_7_length_1600_cov_29.737950  | (2, 234)                        | (2, 234)                       | [(2, 21), (21, 34), (35, 68), (68, 132), (133, 165), (166, 204), (205,   234)] | 74.15                           | 74.15                          | -1         | (199, 1530)                   | (199, 1530)                  | [(1473, 1530), (1313, 1353), (1127, 1228), (770, 963), (594, 693), (395,   513), (199, 284)] | [(1473, 1530), (1313, 1353), (1127, 1228), (770, 963), (594, 693), (395,   513), (199, 284)] | N/A                   |
|                                       | VLNB-5974 | 234          | NODE_10_length_1012_cov_73.506215 | (0, 68)                         | (4, 68)                        | [(0, 21), (21, 34), (35, 68)]                                                  | 81.16                           | 86.15                          | 1          | (159, 546)                    | (171, 546)                   | [(159, 222), (303, 343), (445, 546)]                                                         | [(171, 222), (303, 343), (445, 546)]                                                         | N/A                   |

The useful information is the name of the sequence (Column 4) and the hits match range (Column 12). We will use these information to extract the contig that has the hit and trim the contig to the hit only segment.

Have a code to extract the exon from the fasta and output. The code is tested and works fine:

```bash
python exon_hits_collect.py ../hyb_output/knownmix/ ../output/exon_hits/
```

There are some blank sequences in the FASTA. What happened? (Deal with it tomorrow, really affect the automatic process).
> Update: 09/19/2023 Fixed. The issue was caused by the range of selection. The export of exonerate gives a range can start from 0 but in our script we select [start - 1: end] where when we hit a range like (0, 599) it will extract nothing. A simple fix is taken (and you should know how to :)).

The filtered one looks nice. I neede to generate a mixed tree out of it and see the resolution.

---

I want to separate different exons as different genes, then extract the segment that mapped to the exons.

I would like to edit the remove_nonoverlapped.py to allow a parameter about the size selection intensity. Also, I want to rename the proportion of overlapped parameter.

There are two new side scripts modified from original codes. These scripts are uniquely used in extract exons based on the tsv files. Since most of the steps and thought were recorded before, I gonna show the full pipeline here:

```bash
# A pipeline of reconstructing phylogeny trees using exon information from HybPiper output
# Nan Hu, 09/19/2023

# -----
# This is for testing one gene tree and verifying the pipeline.
# Constants. Please edit them to match the data structure of your own directory. Dont include the last "/" in the path.
hyb_parent_dir=../hyb_output
hyb_species_list=namelist.txt # Equals to the directory list of the HybPiper parent directory
output_dir=../output
testing_gene=5974
expected_tree=true_topo.tre
threads=12
iqtree_dir=~/software/iqtree/iqtree-2.2.2.7-Linux/bin
min_overlap_pct=0.8 # For remove_overlapped.py -i
min_overlap_prop=0.1 # For remove_overlapped.py -p


# Collect exon hits from HybPiper folder based on the information from exonerate_stats.tsv
while read line ;do python exon_hits_collect.py ${hyb_parent_dir}/${line}/ ${output_dir}/exon_hits/ ; done < ${hyb_species_list}

# Create a list contains shared hits among all input species
while read line ; do python find_and_list_exons.py ${hyb_parent_dir}/$line/ ${output_dir}/valid_exons/ ;done < ${hyb_species_list}
python gene_intersection_finder.py ${output_dir}/valid_exons/ ${hyb_species_list} ${output_dir}/shared_exons.txt

# Merge all exon hits of the same gene from multiple input
python merge_exons.py ${output_dir}/exon_hits/ ${output_dir}/merged_exons/ ${output_dir}/shared_exons.txt

# Remove nonoverlapped sequences before alignment to reduce the alignment time
mkdir ${output_dir}/${testing_gene}_test_run
python remove_overlapped.py ${output_dir}/merged_exons/${testing_gene}_exons_merged.fasta ${output_dir}/${testing_gene}_test_run/${testing_gene}_reduced.fasta ${output_dir}/${testing_gene}_test_run/${testing_gene}_nonoverlapped.fasta -i ${min_overlap_prop} -t ${threads} -p ${min_overlap_prop}

# MAFFT alignment and trim over-gapped
mafft --preservecase --maxiterate 1000 --localpair --adjustdirection --thread ${threads} ${output_dir}/${testing_gene}_test_run/${testing_gene}_reduced.fasta > ${output_dir}/${testing_gene}_test_run/${testing_gene}_aligned.fasta
trimal -in ${output_dir}/${testing_gene}_test_run/${testing_gene}_aligned.fasta -out ${output_dir}/${testing_gene}_test_run/${testing_gene}_trimmed.fasta -gt 0.5

# A more stringent composition test to remove sequences having different base compositions from others
python composition_test.py ${output_dir}/${testing_gene}_test_run/${testing_gene}_trimmed.fasta ${output_dir}/${testing_gene}_test_run/${testing_gene}_testPass.fasta

# Realign the rest of sequences
mafft --preservecase --maxiterate 1000 --localpair --adjustdirection --thread ${threads} ${output_dir}/${testing_gene}_test_run/${testing_gene}_testPass.fasta > ${output_dir}/${testing_gene}_test_run/${testing_gene}_realigned.fasta

# Tree construction and comparison
${iqtree_dir}/iqtree2 -s ${output_dir}/${testing_gene}_test_run/${testing_gene}_realigned.fasta -m MFP -bb 1000 -redo
python topo_comp_show.py ${output_dir}/${testing_gene}_test_run/${testing_gene}_realigned.fasta.treefile ${expected_tree}

```
There is a version for batch jobs. Check the input file size before running the job.
```bash
# A pipeline of reconstructing phylogeny trees using exon information from HybPiper output
# Nan Hu, 09/19/2023

# -----
# This is for the full pipeline on all genes.
# Constants. Please edit them to match the data structure of your own directory. Dont include the last "/" in the path.
hyb_parent_dir=../hyb_output
hyb_species_list=namelist.txt # Equals to the directory list of the HybPiper parent directory
output_dir=../output_exons
expected_tree=true_topo.tre
threads=128
iqtree_dir=~/software/iqtree/iqtree-2.2.2.7-Linux/bin
min_overlap_pct=0.8 # For remove_overlapped.py -i
min_overlap_prop=0.1 # For remove_overlapped.py -p


# Collect exon hits from HybPiper folder based on the information from exonerate_stats.tsv
while read line ;do python exon_hits_collect.py ${hyb_parent_dir}/${line}/ ${output_dir}/exon_hits/ ; done < ${hyb_species_list}

# Create a list contains shared hits among all input species
while read line ; do python find_and_list_exons.py ${hyb_parent_dir}/$line/ ${output_dir}/valid_exons/ ;done < ${hyb_species_list}
python gene_intersection_finder.py ${output_dir}/valid_exons/ ${hyb_species_list} ${output_dir}/shared_exons.txt

# Merge all exon hits of the same gene from multiple input
python merge_exons.py ${output_dir}/exon_hits/ ${output_dir}/merged_exons/ ${output_dir}/shared_exons.txt

# A batch run for all genes:
mkdir ${output_dir}/phylo_results
while read gene_name
do
	# Remove nonoverlapped sequences before alignment to reduce the alignment time
	python remove_overlapped.py ${output_dir}/merged_exons/${gene_name}_exons_merged.fasta ${output_dir}/phylo_results/${gene_name}_reduced.fasta ${output_dir}/phylo_results/${gene_name}_nonoverlapped.fasta -i ${min_overlap_prop} -t ${threads} -p ${min_overlap_prop}
	
	# MAFFT alignment and trim over-gapped
	mafft --preservecase --maxiterate 1000 --localpair --adjustdirection --thread ${threads} ${output_dir}/phylo_results/${gene_name}_reduced.fasta > ${output_dir}/phylo_results/${gene_name}_aligned.fasta
	trimal -in ${output_dir}/phylo_results/${gene_name}_aligned.fasta -out ${output_dir}/phylo_results/${gene_name}_trimmed.fasta -gt 0.5
	
	# A more stringent composition test to remove sequences having different base compositions from others
	python composition_test.py ${output_dir}/phylo_results/${gene_name}_trimmed.fasta ${output_dir}/phylo_results/${gene_name}_testPass.fasta
	
	# Realign the rest of sequences
	mafft --preservecase --maxiterate 1000 --localpair --adjustdirection --thread ${threads} ${output_dir}/phylo_results/${gene_name}_testPass.fasta > ${output_dir}/phylo_results/${gene_name}_realigned.fasta
	
	# Tree construction and comparison
	${iqtree_dir}/iqtree2 -s ${output_dir}/phylo_results/${gene_name}_realigned.fasta -m MFP -bb 1000 -redo
done < ${output_dir}/shared_exons.txt

# Testing topology
mkdir ${output_dir}/all_trees

cp ${output_dir}/phylo_results/*.treefile ${output_dir}/all_trees

while read gene_name; do python topo_comp.py ${output_dir}/all_trees/${gene_name}_realigned.fasta.treefile ${expected_tree}; done < ${output_dir}/shared_exons.txt
```

Let's see which method would be better for create a reference tree (and select good genes) for species identification.

