# What kind of data can be used for estimating the mixed species number in the sample?

Why do we need to estimate the species mixed in the sample? Becasue we find out that the parameter choice may be related to the species number mixed in the sample. Also, if we could know the number of species in the mix, we can directly make a truncating in the cumulative similarity table to decide the mixed species.

---

Today I will try to output some numbers that may have correlations to the species number in a mixed sample. The first try will be the contig number that assembled by SPAde. This will be the easiest number we could get.

Then, I want to output the exon frequency from each gene. I will group the exons based on the exon tree. The statistics from this will be how many repeated exon called for each potential species in the mix.

---

I would like to check the individual successful rate in the output. I will use `Mean - 0.2 * SD` as the threshold since it performed the best by empirical data. The reason I did this is because I doubt that there are some individuals hard to be identified through the distance method.

I find that some individuals can be constantly failed due to unknown reason. Actually, if some of them were removed, the result will be very good. Not quite sure if I need to remove those individuals at this moment.

I will run the reduced pipeline to the current data. I will offer a list of mixed sample that exclude those hard to identified individuals.

---

Check QIIME about the standards of cutting down individuals/exon filtering criteria (field standard).

---

Figure out one way to give each species a probability in our final estimation based on 100 sampling table. The idea is, for each value, calculate the delta for `Mean + delta * SD`. Use the delta to check the PPV and NPV in the testing datasets. To do this, I need to generate a table for the confusion matrix by fine steps of delta in our testing datasets. Then, I can assign a number to each species a probability of True Positive or True Negative based on the pre-trained results.

Also, I want to check if there are any interactions between certain species in the candidates. After that, I would like to make another mix that have even number for each candidate in the total pool.

Another thing is, it's time to get back to our original mixes. We can try data 2 & 3 for our known and unknown species once we figure out the probability work. I would like to see the predicted power of the current pipeline.

---

Get some deep thinking of the prediction values in the work of "giving each estimate a probabilty". The key to assess the probability of the being true positive at certain delta value `Pr(is True Positive| K > Mean + Delta * SD)` is PPV (= TP / P). However, the probability of the value to be a true negative at certain delta value `Pr(is True Negative | K < Mean + Delta * SD)` is actually the FDR (= 1 - PPV). It is because the certain value is either being correctly predicted as positive or negative thus the probability of the value being identified as negtive is just the 1 - PPV.

Glad that I can figure this out. I will run the code tomorrow to try one sample. Also, another plan for tomorrow is to redo the mix data 2 and 3 to see the power of the current pipeline.

---

I will estimate the data from our previous dataset, give a probability and try to tease apart how many species hidden in the mix.

---

Got a lot of real data to test on. There are several species that I can identify correctly but forr the Monocot, especially the non-grass ones, there is an issue that they will evenly grouped with all of those references plus stout camphor, water lily, and amborella. I wish if I can figure these out by calculate the node probability. Also, the probability is not quite working here I need to tune up the values.

---

What is the possible method to figure out the samples in monocot? I would like to check the gene trees for those monocot groups to see whether some of the species is missing in the tree. Currently, our pipeline is solely based on the phylogeny. If the phylogeny has issue, the distance method would give wrong prediction in some ways.

---

Don't forget to test the result from different mix compositions!

---

I realize that the different reference should use a different delta value to estimate the probability. This will be a problem for me right now since it limits the usage for custom datasets. 

Let's deal with the issue in those monocots clade first! My current thought is, to category the candidates in the pool, we could first estimate the chance that there is/are monocot(s) in the mix. If yes, then we do something special to it. If not, I am still very confident about our current method in eudicots. 

Before we made this ultimate change, I noticed that the support value for the "problematic" species assignment is pretty low. I think we can add a filtering function in the distance calculation and similarity converting code. The logic will be:
1. Read in a tree from the input
2. Reroot the tree by 'amborella', if not exist then 'water_lily', if both not exist, root by midpoint
3. Check the leaf node that have a string "NODE". Save those as unknown species
4. Check each unknown species. Check the node support to its most recent common ancestor. If the support is >0.75, then record what is the sister species in the tree. If the sister is a clade, record the clade in a form (species1, species2, species3, ...). Notice, if the sister species is also an unknown species then directly goto the second recent common ancestor. If the clade contain an unknown species, remove that from the record.
5. If the support is <=0.75, then go up to the second common ancestor, check the support again, if support is > 0.75, record the clade as the format I showed in step 4. Remove any unknown species in the record too.
6. Repeat step 5 until find a node support > 0.75. 
7. Output the result to a table

After tuning the python code, now I can get the clade that have the most confidence. I will run this for all the trees in a sample output then we can figure out the next step - To use the frequency of these clades to decide the target species.

> The 700+ Angiosperm353 version 2 target file and the alignment is ready now. Super excited to it. 

I will try to make some initial work on this dataset. First, let's start output summary tables.
1. Read in fasta files in current folder. Each file is a gene alignment that named by its file name (without .fasta)
2. Summarize each fasta, how many aligned sequences in the fasta. Output into a output table.
3. For each gene name, for example, >Liliales_Melanthiaceae_Xerophyllum_asphodeloides_1KP-4890 transcript_7797-AFLV, The Liliales is the order name, Melanthiaceae is the family name, Xerophyllum_asphodeloides is the species name, 1KP is the data source.
4. Summarize each species how many genes they showed up in the alignment. Also, when output into a table, add the order and family information as well.

I will take some time to sort the order into a evolutionary reasonable list.

Back to the cladology work. I wish there is a statistical trend to see. Also, I am running a test run for HybPiper using the new Angiosperm353v2.0 for our artificial mix data. I use a full run here to recover all of the genes from the assembly.

A full run for over 700 species in the tree is a huge load for HPCC right now. The HybPiper step is fast but the alignment and tree construction is super slow. I manually checked some trees and find out that the prediction is not bad for our current data. However, I need to change some code for the later work in the python script because the current reference naming system is different from previous 70 species panel.

I submitted the data2, gard10, and corn20 with new fresh full genes runs. I want to see if the pipeline workable with tons of data. I know there will be some codes required to be modified to fit the new reference naming systems but I will do it later. The idea now is to run it and compare the power of one species to be identified in the pool.

---

A week of maintenance in HPCC so I have nothing to run but plan ahead for next week work. I currently have the results for 700+ species full run. I will evaluate the result from those runs then figure out a good way to prune the tree from the 700+ species references. Ideally, I wish to have 2 species for each order unless some tiny orders. Also, I want to give each species hit a series of tags, which will be used in final evaluation. By this, we can estimate the potential hits in different hierarchies. The last thing to do is to thin down the genes/exons in the run to reduce the computational time.

---

Checking the results from the full run. Some of them only have partial results due to the runtime issue. However, I think those would be good enough to have a reasonable estimation in order to look at the predictions. I category the result based on the order, the family, and the genus level. The 80% saff + 20% corn showed great result. At the order level, the Poales and the Asparagales are the only two orders that beyond the threshold. At the family level, we still see Poaceae passes the line and the Asparagaceae followed as the 2nd (but p = 0.058). 

---

Evaluating the results from yesterday and figure out what is missing in the key steps. I think it is reasonable to use different distance calculation methods to treat with different hierarchical search. For a general order search, we could use a more spread algorithm (general distance method, do not need to filter the side hits), while in a specific search we would like to use a more stringent one (the one that only pick up the real closest hits). I will modify the distance to matrix python code to accommodate both requirements by adding an additional parameter.

---

I tested several real mixed datasets. The current conclusion is, for order level identification, our pipeline works well. We just need to tune in the distance acceptance z score in dist2Z.py. For a small reference, we should use lower Z-score or directly use the `--use_flag` to lock on only one hits on the tree. It may not work well in a mix that the composition is super skewed (e.g. 95%:5%). I also found some relationship between the cumulative similarity and the mixed proportion but it is very comparative not quantitative.

I want to write an automatic evaluation code to output the prediction on order, family, and genus level. Then, I want to test a series of parameter settings on distance transformation methods to see which one works better in specific situation.

---

I checked the results from the data3. It is marked as a bad data but by current judging criteria it should still recover something. However, when I checked the results from our pipeline, it showed random results with different parameter choice. I will figure out the problem by looking at the trees and the assembly. I was told that some repetitive regions might be assembled by hybpiper causing the issue.

---

I'd like to put some info from my scratch paper to here.

I noticed a pattern in exploring the data: For those mixed samples that have issue in identification (mostly aore missing taxa), they tend to have less reads mapped to the reference Angiosperm353 ( < 1M reads). However, the 3rd dataset from unknown seed mix and the 4th dataset from pooled soil do have very high mapping rates. I want to see what's wrong with these results by retry these data with new pipeline.

In the meantime, I discussed with Dr. Johnson about some phenomenons found in Haley's results. There are some requirements for HybPiper to correctly assemble the read into supercontigs. First, empirically, it requires at least 100k reads mapped to the target in a single sample to recover 250 genes. A BAM file examined less than this mapped reads number will have trouble in SPAdes. Second, the recovered sequences should target most of the genes with relatively high coverage. Even with a lot of reads mapped to target, if the 353 genes are not well covered, it will affect the number of supercontigs assembled. Third, we need to check if the on target sequences are just nonsense repeats. The MyBaits kit would have chance to enrich a bunch of repeated regions from the library. If these repeats assembled by themselves, it would be no information for the downstream analysis.

To deal with the 3 requirements above, I will check the output from 3rd and 4th dataset which were marked as deprecated data. I wish to find some patterns matches the issue above to explain the reason why they cannot be recovered.

---

I am going to record what I thought and what I did for today.

Currently, the mis-identification of species is mainly because that a random seq that not matching up with any species so it will be assigned to an outgroup with shorter distance. To address this issue, we could 1) reduce the effect of the short distance, which currently solved by ultrametrication. 2) remove weird nodes/sequences that seems to be incorrect. This is done by manipulate the Z-score threshold by selecting the distance in the matrix. However, these unexpected sequences still popped up around some of our results.

I also find a bug in the exon_extract.py which reads the exonerate table for splitting the supercontigs into exons. We should sort the range of all hits before we ID the exons. Instead of fixing that, I tested that if we pool all the exons from the same gene to build a tree, this will no longer be a problem. I tested this thoughts and it seemed work well. I hope this method could somehow reduce the chance that one single piece being randomly placed on the tree.

I explored the difference between fasttree and iqtree. Iqtree is much slower but seems to be more accurate in compare to fasttree. I think we could probably fix a model on iqtree instead of using fasttree for some analyses.

Finally, I found that in the 3rd dataset (pooled seeds), the exon hits is limited with strong similarities to some unexpected lineages. I want to see the true sequences by just assemble the reads directly with SPAdes. I will see what arere actually captures in the pool.