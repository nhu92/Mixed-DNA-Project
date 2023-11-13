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