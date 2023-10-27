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

I talked to Dr. Yibing Zhang yesterday for questions about PCA and clustering methods. I got the thought about the proper ways to display the results and using parameters to evaluate the clustering. I am thinking the way to pile up the exon data into a gene level. Since they are genetic distances so can we simply add them all together?

---

No. The idea that add up the genetic distance won't work because not all of them have all the exons of one gene (actually normally very few of them have). I will go back to use the genetic distance by site but I need to think about a nice and neat way to merge the distance from all exons.

Another thing is to evaluate the merged gene PCA results (or the clustering) by the false positives and false negatives. Although we need to eyeball the PCA/K-means results but we also need to have a statistic support to decide which are our candidate gene by controlling the false rate. 

This also reminds me an issue of figuring out the gene level tree. The gene level tree will generate distance that is not comparable across different contigs that assembled by SPAdes since they may not have the same number of exons. Then, the only way is to calculate genetic distance per aligned loci. I will generate a code to run this process.

```bash
# A pipeline for PCA analysis from phylogeny to PCA including data cleaning and merging.
# Nan Hu, 10/23/2023

# -----
# Constants. Please edit them to match the data structure of your own directory. Dont include the last "/" in the path.
tree_dir=../raw
output_dir=../output
color_taxa=./colorNtaxa.csv
sample_name=knownmix

mkdir ${output_dir}/

# A batch run for all genes:
while read gene_name_shorter
do
        ls "${tree_dir}/${gene_name_shorter}"*"treefile" > ./loop.treelist.txt
        i=1
        while read filename
        do
                python matrix.py -t ${filename} -o ./${gene_name_shorter}.${i}.matrix
                ((i=i+1))
        done < loop.treelist.txt

        python combine_distance.py "${gene_name_shorter}*.matrix" ${output_dir}/${gene_name_shorter}.combined.csv ${sample_name}
        python rm_outliers.py ${output_dir}/${gene_name_shorter}.combined.csv ${output_dir}/${gene_name_shorter}.cleaned.csv
        python mean_statistics.py ${output_dir}/${gene_name_shorter}.cleaned.csv ${output_dir}/${gene_name_shorter}.stats.csv
        python pca.py ${output_dir}/${gene_name_shorter}.stats.csv ${color_taxa} ${output_dir}/${gene_name_shorter}.pca.svg

done < ./shared_genes_shorter.txt

rm loop.treelist.txt
```

This is the pipeline for PCA. It will take the exon tree as the input to calculate the genetic distance. The 2 Component PCA plot is used for selecting candidate genes that matching up the expected grouping. The real decision-making code will be the K-means clustering results. I also want to generate a distance from clustering to show the closed related species to go deep into the clade.

---

Another talk with Dr. Yibing Zhang. It is suggested that we can directly use the large matrix for clustering. The good point is that we won't lose a lot of data but it will be hard to display them on the PCA figure. Another way for clustering is to use PCA1+PCA2 (maybe PCA3?) for clustering thus we could display the clustering onto the figure. The way to decide a matrix (gene) is a good candidate or not is based on a series of parameter. We can try the bulk scoring system by giving 1 point when a reference species is clustered correctly. Or, we could try the scoring system considering the false positives and false negatives.

The PCA results looks generally good except the Sapindales and Malpighales which should group with Fabids but in PCA they are all closed to Malvids. I will run a species tree methods to check the topology of the expected species phylogeny.

Another minor thing, the normalization for PCA generally giving me worse results compared to the non-normalized. I will check these and possibly try normalize and non-normalize methods for clustering.

Today's task in sequence:
1. Generate a species tree of 70 species to check the topology of Sapindales and Malpighales.
2. Clustering 70x70 matrices of several samples to give a try.
3. Clustering all matrices and calculating the index for choosing candidates.
4. Try clustering from PCAs (normalized).

While exploring the clustering methods, I might choose the Spectral Clustering from scikit-learn. The comparisons of all the clustering criterion is displayed below. The input will be the matrices (we need to transformation I will mention it below too) and a cluster number.
![Clustering Methods Comparison](https://scikit-learn.org/stable/_images/sphx_glr_plot_cluster_comparison_001.png "Clustering Comparisons")

As I said, we need a ata transformation before we run the clustering. We got a distance matrix but what we need is a similarity matrix. The scikit-learn suggests a way of  transformation:
```python
similarity = np.exp(-beta * distance / distance.std())
```

---

I worked late night to produce a code for evaluating the clustering results. We used TP, TN, FP, and FN based parameters including ARI, AMI, and V-measure to evaluate the clustering for each gene to the real clustering. These measurements are included in the scikit-learn package. The results can be merged into a file for ranking the genes. The results do not quite make sense since they are not predictable. The talking this morning is focused on the improvement of the pipeline. Here are things we may try:
1. Try a better distance method by removing the gaps from the alignment.
2. Use PCoA which focus on a distance measurement to lower the dimension.
3. Check each dimension to see what they really resolve (e.g. PC1 is to separate Monocots and Dicots).
4. Try different clustering methods including supervised and unsupervised models.
5. Use different precise recall parameters to evaluate the cluster.

---

Is it better to select good gene set or using a whole dataset to summarize the results?

---

Lets try PCoA today. It seems that it needs the least effort. I will also output the distance (PCA - 3 dimensional) of the closest 3 targets to the unknown testing loci and check if those are in the same grouping with the clustering results. I think this is a better criterion than then general non-determinative parameters. The messing grouping hurts the prediction much more than an incorrect identification of the known mix species.

I found a method that using the neighbor nearing method to lower the dimension. It refers to [this](https://scikit-learn.org/stable/modules/neighbors.html). I will give it a try in some of our dataset. It is expected to cluster the neighboring data compared to the PCA and K-means.

