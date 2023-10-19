# Using PCA to screen out candidate genes and using clustering methods to decide unknown pieces
Calculating genetic distance based on some evolutionary model and using the pairwise distance as the input of PCA and K-means analysis to estimate the composition of the unknown mixed DNA.

---

The first thing to do is to generate a pairwise genetic distance matrix from a phylogeny tree. Each exon should have its own tree. Then, we will calculate the mean distance to another node and the standard deviation for variance controlling. I need an algorithm to drop the unqualified node (ANOVA seems to be to loose, and not robust enough). Alternatively, is there any way to display the sd in PCA plot?

The PCA result will be colored by taxonomy groups. The unknown species will be colored differently (maybe highlighted?).

Scikit-learn has the package for pca analysis. I will generate a PCA matrix from it.

