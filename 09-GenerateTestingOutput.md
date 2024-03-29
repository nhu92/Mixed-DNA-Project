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

```bash
python categ_PCoA.py secondmix_refPCoA.csv secondmix_refPCoA_family.csv family
dos2unix secondmix_refPCoA_family.csv  
paste -d, secondmix_refPCoA_family.csv secondmix.Zaxis.csv  > secondmix_3dcoordfam.csv
python contour_optima_normalized.py secondmix_3dcoordfam.csv secondmix_contourmapfam.svg 0
```