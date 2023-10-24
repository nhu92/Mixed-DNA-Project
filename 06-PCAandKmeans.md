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

