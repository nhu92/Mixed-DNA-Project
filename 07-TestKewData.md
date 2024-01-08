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

Let me summarize the current workflow of this testing:

For testing mix:
1. Trim raw reads using fastp
2. Run HybPiper
3. Extract the exon based on the exonerate stats
4. MAFFT --addfragments alignment, trimal, and fasttree GTR+Gamma
5. Run genetic distance matrix, data clean
6. Cumulative genetic distance as Z axis on the contour map
7. Find local optima

For reference panel:
1. Trimal, iqtree -MFP
2. Clean up gene trees, merge
3. Mean genetic distance matrix, some data transformation
4. PCoA lower dimensions, for the X-Y panel

I will compile these steps into a shell script file.
```bash
# This is a pipeline for analyzing raw reads from a mixed sample to the identification.
# This pipeline is to treat with the mixed reads. There is another pipeline to generate the reference panel.

# Declare the constants
threads=64
read1=
read2=
mega353=
proj_name=
ref_alignment=11_final_alignments_and_gff3 # The directory of reference alignment

# ---

fastp -i ${read1} -I ${read2} -o ${read1}.trimmed.fastq.gz -O ${read2}.trimmed.fastq.gz -j fastp.json -h fastp.html
mkdir hyb_output
hybpiper assemble -t_dna ${mega353} \
        -r ${read1}.trimmed.fastq.gz ${read2}.trimmed.fastq.gz \
        --prefix ${proj_name} \
        --bwa \
        --cpu ${threads} \
        -o ./hyb_output

ls -d hyb_output/* | grep -Po "\d+" > gene_list.txt
mkdir ./exon_extracted
while read line; do python split_exon_extract.py ./hyb_output/${proj_name} $line ./exon_extracted 0.8; done < gene_list.txt

input_exon=./exon_extracted

mkdir ./phylo_results
# A batch run to generate trees for all genes:
while read gene_name_shorter
do
        ls ${input_exon}/*${gene_name_shorter}*.fasta > ./exon.list.txt
        i=1
        while read exons
        do
                # MAFFT alignment and trim over-gapped
                mafft --preservecase --maxiterate 1000 --localpair --adjustdirection --thread ${threads} --addfragments ${exons} ${ref_alignment}/${gene_name_shorter}.fasta > phylo_results/${gene_name_shorter}_exon_${i}_aligned.fasta
                trimal -in phylo_results/${gene_name_shorter}_exon_${i}_aligned.fasta -out phylo_results/${gene_name_shorter}_exon_${i}_trimmed.fasta -gt 0.5

                # Tree construction and comparison
                fasttree -gtr -gamma -nt phylo_results/${gene_name_shorter}_exon_${i}_trimmed.fasta > phylo_results/${gene_name_shorter}_exon_${i}.tre
                ((i=i+1))
        done < exon.list.txt
done < gene_list.txt

rm exon.list.txt
# Distance Calc (in progress)
mkdir all_trees

cp phylo_results/*.tre ./all_trees

tree_dir=all_trees
while read gene_name_shorter
do
        ls "${tree_dir}/${gene_name_shorter}"*"tre" > ./loop.treelist.txt
        i=1
        while read filename
        do
                python matrix.py -t ${filename} -o ./all_trees/${gene_name_shorter}.${i}.matrix
                ((i=i+1))
        done < loop.treelist.txt
done < ./gene_list.txt

rm loop.treelist.txt

python dist2Z.py all_trees ./${proj_name}.summary_dist.csv ${proj_name}
python group_sum.py ./${proj_name}.summary_dist.csv ./${proj_name}.cumulative_dist.csv
grep -v ${proj_name} ${proj_name}.cumulative_dist.csv | cut -f2 > ${proj_name}.Zaxis.csv

```

The pipeline of preparing the reference panel. Make sure the panel used is the same as measuring target species genetic distances.

```bash
# This is a pipeline for transforming reference alignment to a PCoA panel.
# This pipeline is to treat with the reference panel. There is another pipeline to generate the mixed input distance matrix.


# Constants. Please edit them to match the data structure of your own directory. Don't include the last "/" in the path.
iqtree_dir=~/software/iqtree/iqtree-2.2.2.7-Linux/bin
ref_alignment=11_final_alignments_and_gff3 # The directory of reference alignment
proj_name=

# A batch run for all genes:
ls ${ref_alignment} > ./allgenes.txt
mkdir ref_tree
while read gene
do
        # Tree construction and comparison
        ${iqtree_dir}/iqtree2 -s ${gene} -m MFP -bb 1000 -redo -pre ref_tree/${proj_name}
        nw_ed ref_tree/${gene}.treefile 'i & b<50' o > ref_tree/${gene}.collapsed.tre
done < ./allgenes.txt

cat ref_tree/*.collapsed.tre | sed 's/_R_//g' | sed 's/_\([0-9]\{4\}\)//g' > ${proj_name}_merged.collapsed.tre
python tree2PCoA.py ${proj_name}_merged.collapsed.tre ${proj_name}_refPCoA.csv ${proj_name}_refPCoA.svg

```

Now, I am moving the testing the 2nd Kew dataset - the "Brassicaceae medley".

The Brassicaceae has 6 reference species in the 70 species alignment. I want to test the power of our pipeline to precisely estimate the species within a small clade. Also, I want to see if the method can successfully exclude the outgroup from the mixed data. I will first trim the 70 species alignment into 6 Brassicaceae species + 1 outgroup (lavender_scallops should be a good one, not too far but quite distinct). Also, I will mix 3 species from Kew tree of life. 2 from different clades of Brassicaceae and 1 from outgroup. 

For the reference trim, here is the species left and their expected phylogeny:
1. arabidopsis
2. rockcress
3. shepherds_purse
4. saltwater_cress
5. schrenkiella_parvula
6. cabbage
7. lavender_scallops

(7, ((1, (2, 3)), (4, (5, 6))));

This topology is also checked from the astral tree generated from the alignment. I will choose two arabidopsis species and one Lamiaceae species. The goal is to test the ability of sort the "best hits" onto arabidopsis and the outgroup, but not on the other species.

---

I end up with using only one species from the Kew - The SRR22519327: Didymophysa aucheri (to saltwater_cress). The rest two reads I selected were from the 1st dataset we had. The Thelypodium (to cabbage) and Salvia (to lavender_scallops).

---

It's the holiday season and seems everyone is slowing down there work! 

I will check the output of the Friday's batch jobs. This will give me more confident in doing the large testing stuff. If this works, I will run another test from identifying the order, towards identifying he family, and end up with smaller taxanomy groups. If this works, it will be super great for me at this moment.

Fixing some pipeline issues. The reference panel has duplicated nodes I am trying to figure out.

Done! I should use mean_distance_genetree.py and the dist2PCoA.py to transform the trees to the PCoA results. I used PC1-PC2 panel since the PC3 is noninformative.

Fixed a folder name extraction command in target distance generator. The current one I used is:
```bash
find hyb_output/${proj_name} -maxdepth 1 -type d -exec basename {} \; | grep -Po "\d+" > gene_list.txt
```

I resubmitted the job. Hope everything comes out fine!

---

Come by the office to check the job status. If the results get fluently generated, I will start another one for Poaceae species. 

I found one of the code is super empty. I fixed it. Now, the pipeline will continue running. This time since it only have 7 references and possible 3 mixed species, it is expected to run shorter time than the 70 species mixed with 8 unknowns.

Another dataset candidate:
Astrebla lappacea (to finger_millet)
Cyrtococcum oxyphyllum (to switchgrass)

The chosen family is Poaceae, with a Poales outgroup Bromeliaceae Ananas comosus (pineapple).
finger_millet
barley
oropetium_thomaeum
rice
bread_wheat
stiff_brome
chinese_silver_grass
switchgrass
biscuit_grass
foxtail_millet
sorghum
maize
halls_panicgrass
pineapple

---

Checking the results. The pipeline should be in a good shape. It is possible to run larger test about the data. Also, I am expecting the new dataset which claimed to have over 700 species well-aligned. By having this dataset, we can easily subset the reference panel to do precisely matching. Let's start with check the result from the run before the holiday.

The pipeline showed some patterns that is both expected or not. The expected species showed up in the final optima but also brought up the unexpected answer. I am thinking about the over filtering possibilities and low representative references panel could be the reason that the reference panel falls apart. I will use the unfiltered matrix to calculate the Z-axis. Additionally, the PC1-PC2 is good to differentiate the in-out groups. The PC2-PC3 should be good to differentiate the in groups down to species.

The Brassicaceae panel only has 6 species, which might be too low to find the correct species. I have hope for the Poaceae panel which has 13 species.

The Poaceae gives a nice prediction overall. The two targeted species are nearly called from the results. The final result is a little bit off but I think it is due to the big polytomy in the Poaceae tree. 

Now, I will search for another Kew dataset for the test in Poaceae. The goal is to estimate the accuracy and the precision. Also, I will trim the 70 species into a simpler one to sort the species in the very beginning. The full pipeline will be submitted tomorrow.

---

Today's tasks:
1. Trim the 70 species into an initial sorting list to quickly sort species into order level
2. Make a checklist for the pipeline. Add the PC2-PC3 command lines to the pipeline.
3. Run the coarse sorting for 2 databases we have (Poaceae and Brassicaceae)
4. Run another Poaceae run with better species identifications.

I finished the pipeline and submitted to the github. Now it should be workable through Github. The guild is [here](https://github.com/gudusanjiao/Mixed-DNA-Project/blob/main/pipeline/contour_map.md). The pipeline works pretty fine currently. It takes less than 2 hours when queuing finished. The species identification on the PC1-PC2 is fairly good to sort the species into their correct places. Since it is very easy to run the command line, I will test more datasets (Why not? :) )

Two more for the selected mix (6th mix):
Oryza australiensis
Sorghum macrospermum

And for the 7th run, the Brassicaceae will be tested in the initial sorting panel.

The 8th run is 2ndMix to eudicots panel. The 9th run is the 3ndMix to trimmed panel v2. The 10th run is the soil(4th Mix) to the trimed panel v2. The trimmed panel v2 is the updated panel to the v1 where the representative species is selected to ensure at least have a reasonable numbers of sp

---

The run is very smooth now. What I can do now is to test the pipeline parameters with multiple dataset, and, with one dataset but fine tuning the knobs. I will use a matrix to test 2 major inputs for one dataset. I will test 4 levels of unknown species + one extreme mix (2, 3, 4, 5, 20+). And I will use different thresholds for contour selection (0, -0.2, -0.33, -0.5). It will be a 5x4 table showing the changes. Also, I will also test the way to spread the reference panel better. I find that if a reference is the odd one on its side, it will only depend on the threshold settings. At least we can spread the panel and increase the searching radius (now is max/20). 

---

The testing plan will be run in several parallel lines. The main line is to test the initial species sorting panel. The species choice is from various clades to testing the sorting process. I will evaluate the True Positives (TP), True Negatives (TN), False Positives (FP), and the False Negatives (FN). The other line will be the ability to run deep lineages searching. I will evaluate the parameter space in Poaceae to see which combination performs better in this situation.

After I searched for the parameters in each lines, I will choose a value from the confusion matrix to calculate the "effectiveness" of the parameter choices. 

The testing sample will be:
1. Plantago, Descurainia
2. Plantago, Descurainia, Astragalus
3. Plantago, Descurainia, Astragalus, Lotus1
4. Plantago, Descurainia, Astragalus, Lotus1, Salvia
5. Plantago, Descurainia, Astragalus, Lotus1, Salvia, ++ Lotus2, Silene, Thelypodium
6. All above ++ Magnoliaceae Liriodendron chinense, (stout_camphor), Zingiberaceae Aulotandra kmaerunensis, (banana), Cucurbitaceae Cucurbitella asperata(cucumber), Malvaceae Dicarpidium monoicum (cotton), Polygonaceae Rumex spinosus(spinach? some Caryophyllales species), Solanaceae Grammosolen dixonii (tomato), Oliveria decumbens (carrot)

Let me submit these jobs on the queue.

All submitted!
```
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
          12418809    nocona 4spref_p      nhu  R      13:46      1 cpu-23-10
          12418780    nocona 3sptarge      nhu  R      16:53      1 cpu-24-16
          12418755    nocona 2sptarge      nhu  R      19:34      1 cpu-26-49
          12418757    nocona 2spref_p      nhu  R      19:34      1 cpu-24-34
          12418182    nocona INTERACT      nhu  R      55:38      1 cpu-23-26
          12418897    nocona 8spref_p      nhu  R       4:37      1 cpu-25-22
          12418843    nocona 5sptarge      nhu  R      10:10      1 cpu-24-3
          12418844    nocona 5spref_p      nhu  R      10:10      1 cpu-25-23
          12418825    nocona 3spref_p      nhu  R      12:04      1 cpu-23-10
          12418824    nocona 4sptarge      nhu  R      12:13      1 cpu-24-2
          12418954    nocona 8sptarge      nhu  R       0:02      1 cpu-24-21
          12418950    nocona 15spref_      nhu  R       0:20      1 cpu-24-20
          12418948    nocona 15sptarg      nhu  R       0:22      1 cpu-24-20
```

Try to submit the second testing line for Poaceae.
1. Oryza australiensis, Sorghum macrospermum
2. Oryza australiensis, Sorghum macrospermum, Astrebla lappacea (to finger_millet)
3. Oryza australiensis, Sorghum macrospermum, Astrebla lappacea (to finger_millet), Cyrtococcum oxyphyllum (to switchgrass)
4. Oryza australiensis, Sorghum macrospermum, Astrebla lappacea (to finger_millet), Cyrtococcum oxyphyllum (to switchgrass), Poaceae: Gigantochloa atter ERR7621555 (to stiff_brome)

Done! Waiting for the results!

---

Generally the run is fine. There is a small issue in the grouping code that remove all the digitals from the species ID, which make the removal of targeted species in the reference panel failed. I will fix this either by changing the code or just simply change the run name of each into a non-digit job name.

I change the code to use "NODE" as the keyword for removing the target hit from the Z-axis file. I reviewed the results then found out that the hazelnut showed higher hits compare to some minor true hits while most of the trees have no closed target sequences associated with hazelnuts. The issue here is because the branch length of the hazelnut is shorter than other lineages. Thus, when there is only one hit to the reference panel, the hazelnut will gain less negative similarities compare to other lineages then been brought up in the final cumulative similarity matrix. To fix this, I came up with a hard modification way to even up the non-related references. I change the `dist2Z.py` code to calculate the `(col_mean - col_sd)` for each target column in the distance matrices from the exons tree. All of the target-reference distance is larger than the value above will be set to 1 (then in similarity transformation they will be transformed into 0.5). This will set all the irrelevant reference into the same number then even out the final cumulative results.

I also changed the final contour map code. The local optima search function used to have a too simple normalizing equation and now was changed to `(z - z.mean()) / z.std()` to perform a true normalizing. I will keep testing this part.

It seems that 0.5 as the threshold for `(z_normalized.mean() + threshold * z_normalized.std())` is good for 2 species in the mix.

Here I update the code need to be replaced for all the testing runs:
```
group_sum.py
dist2Z.py
contour_optima_normalized.py
01-target.sh
b01-rescue.sh
```

I found that the transform the final Z-axis into 0-1 frame then normalize it would give a better operation to considering teh true hits. Also, it seems that the threshold choice is highly related to how many species in the mix. Thus, it is necessary to figure out some methods to estimate how many species in the mixed samples to form a regression to the threshold choice.

---

From yesterday's end day results, the carrot constantly showed up at the boundary of True Positive/False Positive, which caused by the super negative value after normalization of the non-hit references while missing reference carrot were set to 0 which is not quite negative. To deal with this issue, I need to fix the distance to similarity transformation to set the missing data as 0 and the non-hit closed to 0 too.

I temperately removed the normalization in the `dist2Z.py` code to set the missing data and the non-hit reference similar (non-hit will be 1/1000 and missing will be 0 by default). Hope this will solve the high missing data noise.

I ran all the runs again. I find that the current method is nice when the mixed samples is below 5. It will lose some unknown species when increase the species in the mix. Also, it may increase the false positive when the target reference is missing in some exons. I will continue modifying the `dist2Z.py` by removing all the targeted sequences from the reference panel when calculating the column means and standard deviation. To reduce the false positive rate, I will make the judgement statement more stringent (from Mean + SD to Mean + 1.96SD). I will test if this is better enough to solve the issue. If not, I will then change the similarity equation to a let the similarity go weaker.

---

Initial sorting run can be done by bulk similarity + harsh filtering since the reference panel are distinct enough. However, the Poaceae panel is quite similar to each other thus we still need the normalization plus a looser standard to filter the non-hit species.

---

I will run multiple runs for the initial sorting process. The goal is to find a suitable parameter in each stage and output some statistics for checking the accuracy of the pipeline.

---

There will be 4 tools to be developed. 1) A tool to generate a random mix from an input data set; 2) A tool to initializing multiple testing array; 3) A tool to read the output of each run and perform multiple testing on the results; 4) A tool to summarize all the output and plot the result (in R?).

The first tool is done. `mix_generator.py` will read a folder with pre-named FASTQ reads and generate random mixes out of it. The mean mixed number of species is set to be 5 and followed with a Poisson distribution.

I am searching for candidate species on Kew and make a mix out of them. There is a table for the detailed information of these reads and will be used for the final test.

Made 30 candidates from the Kew. Then, I ran this: ```for i in {1..100}; do python mix_generator.py ../reads_pool/ ../reads_mix/; done``` to generate the random mix file for test.

The second tool is done by 3 bash scripts. They will create multiple testing folders, modifying submission scripts, and submit them batch by batch. To reduce the running time, I reduce the 353 target genes to ~80, selected from the most shared genes among 70 species reference panel.

A small issue was found on the phylo tree construction code. It will search for any combination of numbers in the proj_name. I will fix it later but currently it just affect the efficiency.

---

The small issue is fixed. Also, for the initial sorting panel. The mean - 1.96 * sd is necessary. 

I am initializing a code to automatically test different parameters and summarize them into a single table. Then, we can merge all tables and play with some statistics and plotting.

The evaluation code will output a series of values in t he confusion table, including TPR, TNR, FPR, FNR, PPV, NPV, ACC, etc. These values will be used to select a nice delta value for function mean + delta * sd. I will plot some of these values, also identify the bad samples, classify the valid input. Later I will also link these information to the input reads stats to see the coverage.

---

The accuracy of this pipeline is relatively fine. I will calculate a statistic value showing the probability of a true positive to be missed and the probability of a true negative to be called.

I want to check two alternative changes to this pipeline. 1) The way we compare our target species to the reference. We need a matrix that represent a species relationship to others. By comparing this matrix to the reference species matrices we could find a closed one. I will think how differ this is from our original clustering method. 2) Inspired from meeting with Dr. Chris Jackson, what if I filtered the exons that is shorter than 40bp? I do not think the shorter exons should have the same (or even higher?) weights compared to the long exons. Even this may not make the pipeline more accurate, it will at least reduce the computational time.


---

New year! Back to work. I will think about the possibility of the first thought and test the second. (I got my computer reset and I need to reinstall everything to my new system.)