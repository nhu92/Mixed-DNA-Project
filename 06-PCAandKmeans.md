# Using PCA to screen out candidate genes and using clustering methods to decide unknown pieces
Calculating genetic distance based on some evolutionary model and using the pairwise distance as the input of PCA and K-means analysis to estimate the composition of the unknown mixed DNA.

---

The first thing to do is to generate a pairwise genetic distance matrix from a phylogeny tree. Each exon should have its own tree. Then, we will calculate the mean distance to another node and the standard deviation for variance controlling. I need an algorithm to drop the unqualified node (ANOVA seems to be to loose, and not robust enough). Alternatively, is there any way to display the sd in PCA plot?

The PCA result will be colored by taxonomy groups. The unknown species will be colored differently (maybe highlighted?).

Scikit-learn has the package for pca analysis. I will generate a PCA matrix from it.

Made a distance calculating script and a pca script. Also modified the outlier test code for genetic distance matrix. I also found mistakes in taxa assignment in Chris's table. Here is the code for run one sample:

```bash
python matrix.py -t 4471_exon_1_trimmed.fasta.treefile  -o 4471_1.csv
python rm_outliers.py --input 4471_1.csv --output 4471_1.cleaned.csv
python pca.py --matrix 4471_1.cleaned.csv --table colorNtaxa.csv --output 4471_1_pca.svg
```

The result generally looks fine. However, there are very weird malvids that seems to be closed to fabids. I want to check those species later.

---

