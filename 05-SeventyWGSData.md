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

Another thing we discussed is about to assign the correct topology from the tree. In ete3 package, there is a method to claim a node as "Rosids" or "Asterids" etc. Thus, we could test whether our target unknown sequence's relative position.