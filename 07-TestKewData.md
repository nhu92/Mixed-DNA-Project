# Using multiple samples from the Kew database to test the feasibility of our pipeline
3 different samples will be used for testing our pipeline

---

The first dataset I want to sample from the Kew is a widely distributed cluster containing 7 different species from different major clades. The purpose of this work is to test the ability of our pipeline to correctly sort different species into their correct clades. 

The species I chose from the Kew are (what they should be closed to is in the parentheses):
1. Magnoliaceae Liriodendron chinense   (stout_camphor)
2. Zingiberaceae Aulotandra kmaerunensis    (banana)
3. Cucurbitaceae Cucurbitella asperata  (cucumber)
4. Malvaceae Dicarpidium monoicum   (cotton)
5. Polygonaceae Rumex spinosus  (spinach? some Caryophyllales species)
6. Solanaceae Grammosolen dixonii   (tomato)
7. Oliveria decumbens   (carrot)

I am running the HybPiper on the merged reads of these 7 species. After it finished, I will extract the exons based on the exon_stats file as what I did before. Then, I will run the tree method till we calculate the genetic distance. In the meantime, I still need to think about a general way to output a report from the 353 genetic distances matrix to describe the distribution of the predicted species in the mix.

I will use FastTree instead of iqtree to speed up the process. The evolutionary model will be fixed on GTR+Gamma to shorten the model finder processing.

```bash
# A basic usage of FastTree, installed by conda
fasttree -gtr -gamma -nt alignment.fasta > output.tre
```

Also recording the command lines I used from HybPiper output to the exon extractions:
```bash
ls -d */ | grep -Po "\d+" > ../../scripts/gene_list.txt
cd ../../scripts/gene_list.txt
while read line; do python split_exon_extract.py ../hyb_output/kewmix $line ../output_exon_extracted/ 0.5; done < gene_list.txt
```

I found something wrong when we combine all the distance into a single matrix. The PCA combine_distance scripts are problematic too.

---

Currently, I was running the distance on each exons only. The statistics were pulled from the exon level. I will rank the closed related species from each output then merge together to give a probability based on the frequency (and maybe the read size?). 


The output is giving somehow prediction with some variation (it should be variation not mis-identification since no species were found in the wrong related species group). The frequency is not quite related to reads size since tomato should have lowest frequency but on the contrary it has the highest. I would like to think about a strategy to evaluate the variance in the distance matrix. I will try the normalizing way to transform the matrix (to similarity then normalize). Then, giving a probability distribution of a closed related species, what is the probability the unknown being identified as the same species and what is the CI?

---

Yesterday I had a though about how to use the current knowledge to infer the species distance. First, we will generate a 1-D/2-D/3-D panel of reference species. I will use the gene trees distance matrices to form a PCoA panel to locate the 70 species. Then, the genetic distance calculated from each exon (why not gene? because not all sequences have all the exons and the SPAde assembly is not guaranteed to assemble the same species.) to add up to the panel (for example, the 2-D panel is the base of the X-Y and the distance matrix will be added to Z-axis) to make a contour map. All of the local optimals are the predicted species in theory. This method will consider the variance across the phylogeny tree by adding the closed related species' distance but not suggesting a peak. It should be better comparing to just counting the frequency or ranking the best hits since that was affected by the sample strategy and data source drastically. The local optimal will be less affected by the surrounding high calls because we only care about the local peaks.

The first key step for this thought is to find a good way to make the reference panel. The criterion is that we need the panel to be spread enough but still the closed related species should be relatively closer to each other than the distant species. The PCoA should be valid to work on this. I explored the PC1-PC2, PC2-PC3, and PC1-PC3. It seems the last one did the work better but still some of the species are too closed to each other. I want to use some distribution to transform the data before/after the PCoA to make sure the data looks nice. At this moment, it might suggest us to use 3-D base panel but it is very hard to display the result. Thus, I will first explore the 2-D solutions unless they all fail to resolve the expected species.

It looks like a sigmoid transformation with a standardization is nice for our data transformation. From PCoA analysis, the PC1-PC3 panel gave a good distribution of the species which would be useful when I overlay the distance of unknown exons onto it. Now, lets move to generate those unknown exons part.

Generally, I want to apply a normal transformation to each column of unknown exons to make them add-able. To aggregate all the matrices, I need to fix the rows and columns to match up the correct names for the sum up. The output will be used as the Z-axis of the contour map. I will test the normal transformation approach. If it does not work, I might try some other transformation to expand the data differences. 

The normalized matrix should be positive (? maybe- cause it does not affect the contour plot). (The code is very specified need to make it for common runs.)

---

