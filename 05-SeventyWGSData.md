# Using well developed genome sequencing to build a reference gene list
When we have the gene sequence in a good quality, it will be very easy to use mafft --add to align any unknown sequences within a reasonable computing time.

---

Let me generate a species tree of these 70 species? Manually? I need a taxonomy table for them first. 

I checked gene 4471 (the first one). It seems that all of the exon extractions have some hits with the established alignment from the 70 species 4471. However, for long sequences, not the entire exon extracts can be mapped to the alignment. Some fragments performs better than others. The tree is fine right now with tested species (Salvia) grouped with monkey flower (Mimulus ringens) for most of times. This may indicate a distance method would be possible to identify the species.

My current plan is, from exonerate stats file, extract each hit and make phylogeny independently to assess the phylogeny positions. For each tree, calculate the distance of the targeted species to other species and list them out.

Here is a basic thought about this code:
```python
# extract ranges from exonerate stat file, column 6
    ## recognize the structure [(a,b),(c,d),(e,f)]
# check overlaps, then name each range
    ## 90% exon range overlap will be considered as the same exon
# pull the exons from contig file and named by the exon info, column 13
    ## recognize the structure [(a,b),(c,d),(e,f)]
# output the exon with the exon orders, this output will be the input file for mafft --add
```

Write a Python script to do the following work:
```
1. read in a input directory, a gene name, an output dir from argparse. define the basename of the input directory as dataname.

2. extract exon ranges from an exonerate_stats.tsv (under input directory/gene name/data name/, has a header) from the  7th column. The data in this cell looks like [(a,b),(c,d),(e,f)], and there will be multiple extracted exons ranges (a,b), (c,d), and (e,f). Discard rows include and below the row that the first cell is "Hits with subsumed hits removed". print how many exon ranges found for this row.

3. Give the name to each exon range by order with the format "dataname_gene name_exon_number". print exon range names found for each row. there will be multiple exon names in the same row refer to the original rows in tsv file.

4. Check overlaps of these exon ranges, if two exon range have overlap over given percentage (an input), then they will be named with the same name. print the modified exon range names for each row.

5. Using Bio package to extract contigs from input directory/gene name/gene_name_contigs.fasta based on ranges in column 14. The sequence name in the fasta refers to 4th column. The extracted sequence is named by >dataname_exon_name_number. The exon name refers to the new generated column from previous steps. If the 10th column is "1" then the exon names will be assigned to the extracted contigs in the same order. If the 10th column is "-1" then the exon names will be assigned to the extracted contigs reversely.

# Output all sequences from the same exon into a file: genename_exon_name.fasta
# All input and output are handled by argparse.
```

It is really hard to generate the code for this work. The python script named split_exon_extract.py is now working fine. The only concern here is the exon overlapping algorithm will prefer longer exons and will rely on the first exon it hits. However, it won't be a big deal so far. I will continue using `mafft --addfragment` to quickly align exons to 70 species with aligned genes.

---

The program is estimated to run for 36 hours. However, we could use the already generated data to estimate the potential taxonomy relationships of our known mix data. The plan is to take the all the tips named with sample tags and calculate the pairwise distances from the tree. Then, pool the samples from the same node (should be the same assembly from SPAde) of all the exons for each gene to generate an table of average and standard error. Based on the pairwise distances, if one node is not significantly distant from any other node (the outlier), then we will drop this node (or just the exon of this node?).

Working on using some ANOVA to decide either drop an exon or collapse multiple exons from the same species.

Dr. Johnson introduced a way to shorten the executing time. We could concatenate the exons from the same contig and assume they are from the same species (a little bit risky here but can be solved by later processes). Then, we only need to run the tree for each gene not for each exons. Also, we could run MFP of `iqtree` to find a best evolutionary model for each gene. By doing this, we could save time on running MFP for each exons and will also be consistent when we calculate the genetic distances.

Another thing we discussed is about to assign the correct topology from the tree. In ete3 package, there is a method to claim a node as "Rosids" or "Asterids" etc. Thus, we could test whether our target unknown sequence's relative position. It is the ete3 package [node annotation](http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#node-annotation) feature.

---

I have a general thought about this work now. I will extract the distance from each exon of each gene. Then, merge the distance information from the same Node (contigs assembled by SPAde) by mean and standard deviation. Also, I will use the exon tree from each gene to generate an Astral species tree. Then, the two input file (astral tree & distance statistics) is going to be used in a python script that plotting two information (tree on the left side and a heatmap for the statistics). Additionally, the exon will be filtered based on the ANOVA. I will also summarize the highest hit(s) (top 10?) among all genes to give an estimation of which species it could be in the mixed sample.

---

Today's task:
1. Refine heatmap code. The color gradient need to be changed. 
2. Remove outliers from all exons, based on the ANOVA statistics, should be very loose.
3. Pool the exon results from each node, calculate the mean and sd, format them into the input form. Run heatmap.
4. Collect the outliers from each node, summarizing the species as the candidates for the supportive candidates.

I really need 128 CPUs to run the exon splitting job. I hope it will reduce the runtime into 24 hours. The previous job didn't finish so I could only use around 900 trees to test the result. However, it is not bad at all since it will still have around 150 genes.

Some codes to clean up the tree and form a astral ref as input tree.
```bash
ls | grep -Po "^\d+" | sort | uniq > done_list.txt
while read line ; do cat ../${line}* >${line}_merged.tre ; done < ../done_list.txt
while read line ; do sed -i 's/_R_//g' ${line}_merged.tre; done < ../done_list.txt
while read line ; do python clean_reroot.py --tree ${line}_merged.tre --start_str knownmix --output ${line}_trimmed.tre ; done < ../done_list.txt
while read line ; do astral -i ${line}_trimmed.tre -o ${line}_astral.tre ; done < ../done_list.txt


ls ../*.tre > tree_list.txt
# Node distance calc
while read line ; do python ../../../script/distance_calc.py -t $line -n NODE -o $line.csv ; done < tree_list.txt
mv ../*.csv ./

```

---

The astral tree thought is not valid. The suggested way to generate the tree and the distance is directly concatenating all of the exons. With this way, we do not need to sort the exon and generate exon based table for distance statistics. This approach will be much less computational intense since we only have ~300 genes need to run. The plan today is to generate a concatenated exon with 10 gaps in between. 

1. Generate a code to extract exon and concatenate them together with 10 "-".
2. Try mafft --add to check the alignment for multiple non-overlapping sequences.
3. Run iqtree -m MFP for each gene of 70 species to figure out the best model - in order to skip model selection step for each gene tree ith known mix.

```
1. read in a input directory, a gene name, an output dir from argparse. define the basename of the input directory as dataname.

2. extract exon ranges from an exonerate_stats.tsv (under input directory/gene name/data name/, has a header) from the  7th column. The data in this cell looks like [(a,b),(c,d),(e,f)], and there will be multiple extracted exons ranges (a,b), (c,d), and (e,f). Discard rows include and below the row that the first cell is "Hits with subsumed hits removed". print how many exon ranges found for this row.

3. Using Bio package to extract contigs from input directory/gene name/gene_name_contigs.fasta based on ranges in column 14. The sequence name in the fasta refers to 4th column. Concatenate the extracted sequences from the same row into an output, the gap between two interval ranges are filled by 10 "-"s. The extracted sequence is named by >dataname_{gene name}_{sequence name}. 

# Output all sequences from the same exon into a file: genename_exon_name.fasta
# All input and output are handled by argparse.
```

The `mafft --addfragment` could not correctly align the concatenated exons properly as the separate exons. I had to switch back to original exon way. However, I removed files that is too large for `mafft` and `iqtree` to save some time on the queue.

---

I will finish the statistic code today for the per exon pipeline. The remaining stuff are: 1) drop outliers code; 2) statistics table transformation; 3) standardize the mean; 4) refine the plotting code; 5) cross gene summarization code.

---

New work on chloroplast genome!

---

I am exploring the possibility of use the full length of exons. I found the alignment sometimes is good and sometimes is bad. The tree performance is weird with a long branch length byt sometimes the position is fine. The already aligned data is intron-free so I should not add "-" or "N" between exons when concatenating them. However, I found that the exon information in 70 species is different from the exon identified from exonerate. I would like to explore a little bit more to see what happening there.

---

It's the Herbarium Day. I found an interesting paper just published! [Ancient DNA is preserved in fish fossils from tropical lake sediments](https://doi.org/10.1111/mec.17159).

---

Continuing working on the coding. The current stage is to sort the statistics and output plots. I will think about how to make a better summarization and additionally refine the way of alignment. The exon alignment works fine but it takes too long to run. Either we could align the best model for each gene to save time for model finder process. Or we could try concatenate the exon and try the alignment. I found that exons called by exonerate are shorter than the well-aligned, which is weird to me. I will spend some time to explore it.

---

Yesterday I only finished the combining part. It is harder than I thought. I will finish the statistics half today because the output remover code is already good enough (through ANOVA test alpha = 0.001). Also, I will get the exons extracted from 8 known species (data 1) to examine the distance result by mix the exons to the alignment ith 70 species and known mix. I will use the `intronerate.gff` as the exon information to extract exons from the `{gene}_intronerate_supercontig_individual_contig_hits.fasta`. I checked one gene and the exon that is numbered as the exonerate output and the exon extraction code. And I manually aligned the exons that extracted from `intronerate.gff` to the exons from the knownmix and they aligned pretty well to each other so I think this is workable. I will try to finish the statistics code and then turn to it.

Things went nice today. I did not spend a lot of try runs for the statistics code. I could run the statistics pipeline through all the genes to generate tons of figures! Also, there is a change for a PCI analysis by using the genetic distance (may need transformation since the closed related species have smaller number) from all reference species. 

This command line is to run the exon extraction code:
```bash
 while read gene; do while read species; do python known_exon_extract.py ../../../202309/mixedDNA/1st_8indv/hyb_output/${species}/${gene}/${species}/intronerate/${gene}_supercontig_without_Ns.fasta ../../../202309/mixedDNA/1st_8indv/hyb_output/${species}/${gene}/${species}/intronerate/intronerate.gff ${gene} ${species} ../output/ ; done < ../../../202309/mixedDNA/1st_8indv/script/namelist.txt ; done < ../../../202309/mixedDNA/1st_8indv/output_exons/shared_exons.txt
 ```

 ```bash
 # A typical sample run
for i in $(find *treefile); do sed -i 's/\.exon[0-9]*//g' $i; done
cat *treefile > 4471.gene.tre
sed -i 's/_R_//g' 4471.gene.tre
python ../../../70species/output_exon_split/partial_trees/merged_tree/clean_reroot.py --tree 4471.gene.tre --start_str knownmix --output 4471.cleaned.tre
astral -i 4471.cleaned.tre -o 4471.astral.tre

for i in $(find *treefile); do python ../../script/distance_calc.py -t $i -n NODE -o ${i}.csv;  done
python ../../script/combine_distance.py --input-pattern "*treefile.csv" --output-file ./4471.combined.csv --remove-string knownmix
python ../../script/rm_outlier.py --input 4471.combined.csv --output 4471.cleaned.csv
python ../../script/noded_stats.py 4471.cleaned.csv 4471.stats.csv
python ../../script/rm_outlier.py --input 4471.combined.csv --output 4471.cleaned.csv
python ../../script/noded_stats.py 4471.cleaned.csv 4471.stats.csv
sed -i 's/,,/,0,/g' 4471.stats.csv
python ../../script/tree_heatmap.py 4471.stats.csv 4471.astral.tre 4471.heatmap.svg

 ```