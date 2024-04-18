# Let's make some nice plots for the poster.

I will unify the plot generation using R::ggplot2. 

---
Here is the note from the meeting with Dr. Johnson.

Notes for Poster
1. Use of 800 species as flowering plant phylogeny reference
2. Development of Z-score method for assessing similarity across many loci
    a. topo-map figure, example tree(s)
    b. classifier logic (using topology and the MRCA support values)
3. Positive/Negative predictive value analysis (include definitions)
    a. how many reference samples we need to ID to order: number of genes and number/phylogeneitc spread of taxa
    b.mixture numbers versus accuracy - less “busy” figure by just choosing mixes of 3, 6, and 11
4. Specificity/Precision analysis (minimum reads to ID)
5. Runtime analysis

The first figure would be the regeneration of the topo map. This figure was removed from the current analyzing pipeline. I will recall the code and try to add some reasonable edition to it.

This figure requires the astral tree of the testing reference. Since we use the new A353v2, I would like to evaluate the tree before the plot generation. I could reuse the pipeline from 02-reference.sh without generate the contour map. The pipeline should finish before the exon tree generation process.

Find some old code file and pipelines to do this. I may remove the SD for a reason.

---

I will also generate a contour map for the result. I had never used that plot for the 100 species map. I think it might be messy but I will check.

Glad to pick up some old scripts generated before. I almost forget them already. The contour map will be manually generated from the two inputs: the reference PCoA result and the cumulative similarity result. I might need to fix the PCoA result by making the seq name into specific hierarchies. For example, we need the level at orders/families. Also, the final contour map does not like a very long descriptional naming.

This is the small pipeline to generate the expected contour map. We do not need to run the ref_panel pipeline every time for this work. If the reference file does not change (which means, both the number of genes and the reference species), we can use the same `_refPCoA.csv` as the input file.

```bash
proj_name=
level=
python categ_PCoA.py 113_30_refPCoA.csv 113_30_refPCoA_${level}.csv ${level}
dos2unix 113_30_refPCoA_${level}.csv
dos2unix ${proj_name}.Zaxis.csv
paste -d, 113_30_refPCoA_${level}.csv ${proj_name}.Zaxis.csv > ${proj_name}_3dcoord${level}.csv
python contour_optima_normalized.py s${proj_name}_3dcoord${level}.csv ${proj_name}_contourmap${level}.svg 0
```

---

I will start to rerun the 100 tests from the Kew Garden database. Should I remove the individiuals that we found some issues with? Or I can replace those species with other ones from the Kew? Also, I wish to sample each species with even times to avoid bias. Also, I want to generate a super mix with reduced reads number.

I will check each species reads number, then, I will thin down the read number for each species into 150k reads to avoid a huge input file. 

---

I need to generate some example trees for the display. It is better to color them in a good way. I probably need to show the heatmap with the topology on side. 

---

This is the process of generating new sets of 100 mixes.

I will generate 40 3x mixes (120sp, 6/sp), 40 6x mixes (240sp, 12/sp), and 20 10x mixes (200sp, 10/sp). The detailed processes is divided below:

1. Thin down the reads into 1000k (very sufficient) for each input.
2. Tablize input name to make mix, each column is a random redistribution of 20 names
3. Merge read according to the random table.

The code for making the mix. 
```R
namelist <- c("01x", "02x", "03x", "04x", "07x", "08x", "09x", "10x", "12x", "13x", "15x", "16x", "18x", "19x", "22x", "24x", "25x", "26x", "28x", "30x")
mix3 <- t(replicate(40, sample(namelist, 3, replace = FALSE)))
mix6 <- t(replicate(20, sample(namelist, 6, replace = FALSE)))
mix10 <- t(replicate(20, sample(namelist, 10, replace = FALSE)))

```

---

Submitted 100 jobs on 10 genes tests. I expect some of them failed due to limited gene number. I will expand those tests to more genes tomorrow.

---

This is the new evaluation pipeline for the new 100 tests:
```bash
while read line; do python cumu_order.py ../gene10/${line}/${line}.cumulative_dist.csv ../eval_gene10/${line}.cumulative_dist_order.csv; python evaluate.py ../eval_gene10/${line}.cumulative_dist_order.csv candidates_order.tsv ../eval_gene10/${line}.order_eval.csv; done < mix_namelist.txt
```

