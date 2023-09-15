# Any ways to extract useful information from the output?

---

I am back to work. It's a beautiful Friday. I gonna check if everything is ready for a HybPiper run.

Here is the regex I used to generate the namelist.txt. However, some of the species name need to be corrected manually.

```bash
ls | grep -Eo "(^[^.]*)"  | sort | uniq > namelist.txt
```

Job submitted. I created 8 jobs and each of them requested 64 CPUs. Hope our HPCC staff won't be mad at me.

> In the meantime, I was preparing another run of HybPiper using known mixed samples.
```bash
# similar command lines but do not need to run the while loop
#!/bin/bash
#SBATCH -J hyb_asbl_known_mix
#SBATCH -p nocona
#SBATCH -o log/%x.out
#SBATCH -e log/%x.err
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=3G

hybpiper assemble -t_dna ../raw/mega353.fasta \
        -r ../raw/SSBseq002.trimmed.R1.fastq.gz ../raw/SSBseq002.trimmed.R2.fastq.gz \
        --prefix knownmix \
        --bwa \
        --cpu 64 \
        -o ../hyb_output
```

Back from lunch. I am now exploring intronate folder of the hybpiper output, mainly focusing on the .gff file generated from the gene prediction and alignment based on the supercontigs. The supercontigs are the contigs that assembled from SPAdes concatenated together with 10 "N"s in between (under `hyb_output/<PREFIX>/<GeneNum>/<PREFIX>/sequences/FNA/`). As discussed with Dr. Johnson, he infered that the exon boundaries of the targeted genes for most of angiosperms are pretty much conservative, which means we can extract the shared exons from these supercontigs and reconstruct a gene tree for each exon.

The base gene tree will be only the 8 species we know. This tree will be the base for testing the known mix (artificially mix 8 species samples) data whether it could correctly grouping with their species.

I am going to write some python scripts to fetch info from each .gff file within each species gene by gene. This next level work I think worth a separate chapter.

3. How to extract exons/genes from our supercontigs?

The input file will be under `hyb_output/<PREFIX>/<GeneNum>/<PREFIX>/intronerate/`. I will try to extract the common regions from all genes in 8 species.

(It is Monday again. My nose was bleeding and runny at the same time this morning. Not feel so good right now)

However, after discussion with Dr. Johnson last Friday, It is not suitable to use the extronerate output as the input guidance for extracting all the sequences because it only offers one result while a mixed DNA sequencing contains multiple potential DNA from different species. We need to go back to use an aligner (maybe mafft?) to align all of the assembled contigs.

Before the alignment algorithm considered, we need to fix the contig name from the hybpiper. The hybpiper did not assign proper name to each contig it assembled. Most of them are named "Node #". We want to fix it based on the species/input data name and the gene target name plus an order. Here is the format I would like to have: >species_gene#_order.

I gonna write a python script to substitute the sequence name. It will have 3 steps. First, copy all the contig sequences files into a working directory and rename it with species_gene.fasta; Second, create a substitution list, 1st column will be the old names, 2nd column will be the new names; Third, use SeqIO to change fasta files with new sequence names. (Thanks ChatGPT!)

> I prompt this in ChatGPT3: Give me a python script to: generate a list file of folder names from a given directory; copy fasta files named like "foldername_contigs.fasta" to an input directory; substitute the sequence title of each fasta title in the input directory a pattern with "input directory name_folder name_an order number started from 1". You may use some package like SeqIO.

After several rounds of modifying and prompting, the tool works perfectly now. It is available as [change_seq_names.py](https://github.com/gudusanjiao/Mixed-DNA-Project/blob/main/change_seq_names.py). This tool can be used as shown below:
```bash
# an example to apply to astragalus folder where stores the HybPiper outputs from astragalus input reads.
python change_seq_names.py ../hyb_output/astragalus/ ../output/
```

The output will be in the `../output` folder with the naming pattern "gene#_species_contigs.fasta". In each FASTA file, the sequence name has been changed in a pattern of ">gene#_species_order#".

I will use this tool to extract the fasta file and rename them very quick. Remember that we had generated a namelist file for the while loop.
```bash
while read line; do python change_seq_names.py ../hyb_output/$line/ ../output/assembled_contigs ; done < namelist.txt
```

It's time for lunch but I am not hungary right now. I have a feeling that my runny nose is ready for a big comeback. Let me run the loop and take a nap.

Come back from the rest. I feel better now. The AC of the fifth floor is really really strong btw.

4. MAFFT alignment: Find the common region of all species in each gene.

The basic idea for this part is to take the assembled contigs mapping to a target reference. This reference can be the subset of mega353.fasta (which is pretty shattered), the FNA output from HybPiper, or the supercontigs predicted by Exonerate. I will try all of them to see which one could be a better candidate using the 8 species data.

Before we running mafft for alignment, several steps remains to be taken. We need to filter out genes that failed to be assembled in some of the species. We also need to extract our target reference file from the HybPiper (we can use hybpiper retrieve_sequences, so I will try this first).

More Python scripts incoming! I love ChatGPT!

> This is the prompt I used to generate the frame of the code: Give me a python code to: find all the 1st level folder names in the source directory; within each folder, check if the file with the pattern "folder name_contigs.fasta" exists; list folder names that have that file output into a txt file named "{source directory name}_assembled_genes.txt"; need to have an input source directory and an output directory using parse.

Works perfect! I named it [find_and_list_assembled_genes.py](https://github.com/gudusanjiao/Mixed-DNA-Project/blob/main/find_and_list_assembled_genes.py). Here is a sample run:

```bash
while read line; do python find_and_list_assembled_genes.py ../hyb_output/$line/ ../output/valid_genes/ ; done < namelist.txt 
```

Then, I "requested" a code that find the shared gene names across all of the gene lists generated from last step. The code named [gene_intersection_finder.py](https://github.com/gudusanjiao/Mixed-DNA-Project/blob/main/gene_intersection_finder.py)

> This is the prompt I used to summarize the result from gene lists: Give me a python code for: read in a list of files that in a pattern "{species_name}_assembled_genes.txt", the species_name are input from a list through parse; find the shared gene names among all the list; require an input directory of the list files in parse; require an output as "shared_genes.txt" and an output directory in parse.

This script takes 3 input: the directory with gene lists, the name list file, the output directory.

```bash
python gene_intersection_finder.py ../output/valid_genes/ namelist.txt ../output/
```

Now, we have the target genes that we would like to align. I think one more python code would run the alignment pretty well. For this alignment, as suggested by Dr. Johnson, we could use `mafft --add` to speed up the alignment since we do have a lot of sequences to run. 

Retrieving fasta files from HybPiper using `hybpiper retrieve_sequences`:
```bash
hybpiper retrieve_sequences -t_aa ../raw/mega353.fasta dna --sample_names ../script/namelist.txt
```

It's a little bit disappointing that I could not test the alignment part. However, I have all the files that ready for the alignment test, wish me a good luck tomorrow. (Today's D20 is 17, not bad.)

---

Tuesday! Cold one! (D20 = 4, let me wish some luck for today)

Today I want to explore the possibility of using python to run `mafft` for a sequence alignment.

The first step I would like to have is to pool all the candidate sequences together. The retrieved sequences are served as references. The python code is named [merge_ref_contig.py](https://github.com/gudusanjiao/Mixed-DNA-Project/blob/main/merge_ref_contig.py)

```bash
python merge_ref_contig.py fna_dir fasta_dir output_dir gene_list_file
```

> The initial prompt I used: Give me a python code that takes the input list of gene names, search FNA file name that start with the gene name in the first input directory, merge with the fasta file name that start with the same gene name, do it for each gene names in the input list. output as "genename_ref_contig_merged.fasta". all input and output directory are input from argparse.

I successfully merged the reference FNA and the contigs that assembled by SPAde. I would like to try some sample trees using `iqtree`.

> I tried some target output from some genes, tested 3 genes trees does not exactly have the same topology as they should. The expected tree in Newick form should be: ((plantago, salvia), ((((lotus1, lotus2), astragals), silene), (descurainia, thelypodium)))

I also tried trees from some merged data. The alignment was messy (with 0 gap-free sequences), I need to find some methods trim/remove some of the contigs and at least leave one contig per species. There should be a tool to do this or I have to make some criteria to filter the results before running the mafft.

I am using gene 7361 to try the filtering. Currently, I am trying to do pairwise alignments to filter out the contigs that has less than 400bp overlapping with at least 10% of the input contigs for gene 7361. This is a really harsh fliter, with 11/653 sequences passed. The mafft alignment discovered 56 non-gap loci which is not bad for a tree.

The code for pairwise alignments is called [remove_nonoverlapped.py](https://github.com/gudusanjiao/Mixed-DNA-Project/blob/main/remove_nonoverlapped.py). I tried to make it run in multithreads but I am not sure if it is actually working properly. It takes a while for the job to be done.

For testing the filtering, here are the bash commands and software I used:
```bash
python remove_overlapped.py ../output/merged_ref_contigs/7361_ref_contig_merged.fasta 7361_reduced.fasta 7361.nonoverlapped.fasta --min_overlap 400 --num_threads 128 --filter_intensity 0.1
mafft --preservecase --maxiterate 1000 --localpair --adjustdirection --thread 128 7361_reduced.fasta > 7361_aligned.fasta
~/software/trimal/trimal/source/trimal -in 7361_aligned.fasta -out 7361_trimmed.fasta -gt 0.5
~/software/iqtree/iqtree-2.2.2.7-Linux/bin/iqtree2 -s 7361_trimmed.fasta -m MFP -bb 1000 -redo
```

Filtering results are like below, with 653 raw sequences (if too many sequences remained, there will be no alignment and tree construction):
| Min Overlapping | Min Overlapped Sequences Proportion | Removed Seq Counts | TreeExpTopo(Y/N)  |   |
|-----------------|-------------------------------------|--------------------|-------------------|---|
| 400             | 10%                                 | 642                | Y(-1species)      |   |
| 400             | 20%                                 | 644                | Y(-1species)      |   |
| 350             | 10%                                 | 628                | N                 |   |
| 350             | 20%                                 | 635                | N(&-1species)     |   |
| 300             | 10%                                 | 387                | -                 |   |
| 300             | 20%                                 | 553                | -                 |   |
| 300             | 30%                                 | 601                | N(messy)          |   |
| 250             | 20%                                 |                    |                   |   |
| 250             | 30%                                 |                    |                   |   |

I noticed that there is a composition test at the very beginning of the `iqtree` run. I think this is very important to filter our data since a lot of 'N's and gaps are generated through HybPiper and alignment. I would like to remove those sequences that failed to pass the chi-sq test (have a strong deviation from average base composition).

The composition test is to count the composition of the 4 bases in a sequence and perform a chi-square test between the individual sequence composition and the average of all the input aligned sequences. The significant differed one (p<0.05) will be removed from the input. The `iqtree` one consider 'N' as the wild card for calculating p values. However, the script I have here will treat 'N' as gaps thus would be more rigorous than `iqtree`. The python code is named as [composition_test.py](https://github.com/gudusanjiao/Mixed-DNA-Project/blob/main/composition_test.py). There is an example of usage below.

The testing code is updated into:
```bash
python remove_overlapped.py ../output/merged_ref_contigs/7361_ref_contig_merged.fasta 7361_reduced.fasta 7361.nonoverlapped.fasta --min_overlap 350 --num_threads 128 --filter_intensity 0.1
mafft --preservecase --maxiterate 1000 --localpair --adjustdirection --thread 128 7361_reduced.fasta > 7361_aligned.fasta
~/software/trimal/trimal/source/trimal -in 7361_aligned.fasta -out 7361_trimmed.fasta -gt 0.5
python composition_test.py 7361_trimmed.fasta 7361_testPass.fasta
~/software/iqtree/iqtree-2.2.2.7-Linux/bin/iqtree2 -s 7361_testPass.fasta -m MFP -bb 1000 -redo
```
```powershell
# the output from the composition test
Sequence 7361_astragalus_1: Passed, p-value: 0.8928951711975892
Sequence _R_7361_astragalus_2: Passed, p-value: 0.6250132791405318
Sequence 7361_astragalus_3: Failed, p-value: 0.032271787690734445
Sequence _R_7361_astragalus_4: Failed, p-value: 8.807631565729831e-09
Sequence 7361_astragalus_5: Failed, p-value: 0.0052726599343652075
Sequence 7361_plantago_1: Failed, p-value: 9.545921467408154e-06
Sequence 7361_plantago_2: Failed, p-value: 0.001631472660128785
Sequence _R_7361_plantago_3: Failed, p-value: 1.0251301030931275e-09
Sequence _R_7361_descurainia_1: Passed, p-value: 0.27583618458958103
Sequence 7361_lotus.str_1: Passed, p-value: 0.7661602754092356
Sequence _R_7361_lotus.str_2: Failed, p-value: 0.004850103401797057
Sequence 7361_lotus.str_3: Failed, p-value: 0.0011918671284876204
Sequence _R_7361_lotus.sal_1: Passed, p-value: 0.7890423551406938
Sequence 7361_lotus.sal_2: Failed, p-value: 0.010204448451551952
Sequence 7361_thelypodium_1: Passed, p-value: 0.06252181970377735
Sequence 7361_thelypodium_2: Failed, p-value: 0.006089466908471762
Sequence _R_7361_thelypodium_3: Passed, p-value: 0.29135159596677596
Sequence _R_7361_thelypodium_4: Failed, p-value: 3.4131290440812136e-07
Sequence 7361_salvia_1: Failed, p-value: 0.017341859225610646
Sequence _R_7361_salvia_2: Passed, p-value: 0.24217166981498744
Sequence 7361_salvia_3: Failed, p-value: 3.8960814550266134e-10
Sequence 7361_salvia_5: Failed, p-value: 7.944164883052365e-05
Sequence _R_7361_salvia_6: Failed, p-value: 0.0044219420625107335
Sequence 7361_salvia_7: Failed, p-value: 0.00016060252759442475
Sequence _R_7361_salvia_8: Failed, p-value: 1.0213913443131421e-22
```
This gene (7361) totally miss the I would like to run these tests for a larger set of genes. I would like to random choose 10 genes and check their output tree.

```bash
while read line
do
        python remove_overlapped.py ../output/merged_ref_contigs/${line}_ref_contig_merged.fasta ${line}_reduced.fasta ${line}.nonoverlapped.fasta --min_overlap 350 --num_threads 128 --filter_intensity 0.1
        mafft --preservecase --maxiterate 1000 --localpair --adjustdirection --thread 128 ${line}_reduced.fasta > ${line}_aligned.fasta
        ~/software/trimal/trimal/source/trimal -in ${line}_aligned.fasta -out ${line}_trimmed.fasta -gt 0.5
        python composition_test.py ${line}_trimmed.fasta ${line}_testPass.fasta
        ~/software/iqtree/iqtree-2.2.2.7-Linux/bin/iqtree2 -s ${line}_testPass.fasta -m MFP -bb 1000 -redo
done < subsample_genelist.txt
```

I am evaluating the final tree topology. I found some pattern that is wordy to explain. I need to find someone to discuss about it.

---

Wednesday morning after a long meeting with Dr. Johnson. Get some ideas to filtering the merged contigs.

- Try treeshrink on each gene tree to remove bad branches (too long ones).
- Explore `exonerate` output from the `HybPiper` and fetch the contigs that have hits on the exon report.
- If can't find proper outuput, run `ryo` options of exonerate to form a useful & reliable contig list.
- After the composition test, rerun `mafft`, since the one cause the alignment issue may be gone.
- Ideas about how to compare the topology of two trees. a) Fetch the sequences names from the query tree. Make the expected one of those names. Test topology; b) Do the monophyly test; c) Calculate the variance of the pairwised phylogenetic distance.
- To get "reliable" results from the mixed-DNA reads, other than using HybPiper (which sort the reads to each gene then assemble to the target), we could try to do the assembly of all the reads (overlapping method) first, then assign the assembled sequences to each gene (hidden markov modeling), which may have better output for an unknown mixed samples. (But this won't be suitable for a qantitative test right?)

Get a lot to do today! I will first try my approach on the known-mixed samples. Then identify around 5 genes that has a expected topology of their gene trees and extract those sequences and merged with mixed contigs and then output a tree. The rest of the time today I will try treeshrink & exonerate information about matched contigs. (D20 = 20, YAY!)

I ran the test run for 8 species data using 400 as the minimum overlapping. The output trees still have a lot of bad sequences unfiltered. I wonder if I need to modifiy some code to set up the minimum overlapping bases based on the minimum length and the average length of the input file. I have the log file from the test run so I could figure out a nice way to filter the data.

There are some long branches can be trimmed. But that is not the only problem here. Some phylo groups have weird topologies. I would like to explore what is happening there.

I also ran the pipeline for known mixed data. I would like to merge them to a 'reliable' set of sequences and see if any information could be extracted from the topology.

Since mentioned the topology, I have a rough thought about how to compare the phylogeny input to the expected phylogeny. Basically, I would like to form a guide tree (a sure one). Then, based on the tip names from the input tree, let the branch grow. Also, remove all the affixes of the tips (only species left). Finally, find an approach to compare the topology of two trees with the same elements.

---
Thursday morning (D20 = 18). I reviewed trees generated from the last test. It turned out be have 4 problems. The most common one is the outliers having long branches. This is relatively easy to handle since we could be more stringent in filtering or trim the final tree using some packages. The second issue is about the unexpected topology. This was happening mostly about the relative phylo position of Silene. This is chanllenging to explain but easy to deal with just select the gene with the expected topology. The third issue is the low supporting values. Most of the low supports appeared with the extra long branches thus I hope this will be solved by thinning down the input sequences. The last issue might be the most problematic one. There are some mixed structure in the gene tree. Especially happened between two Lotus species. In some cases, it can affect two Brassicaceae species, too. This might randomly affect the resolution of our approach to identify species.

I decide to modify the filtering code to take the input file size as consideration. I am not quite worried about the false negatives right now since we do have a lot of candidate genes. My current workflow would be: 
- Filter data considering the input size of the sequences
- Mafft alignment 
- Trimal at -gt 0.5
- Composition test to remove sequences
- Take the remain sequence names, fetch from the filtered FASTA, redo the mafft
- Redo trimal
- Redo compsition test to examine
- Iqtree
- Treeshrink to remove outliers.

I tried `treeshrink` in both alpha = 0.05 and alpha = 0.1 threshold. The gene tree of 5168 has a long branch '_R_thelypodium_9' but it was not trimmed. I will rely on the data filtering side more at this moment.

```bash
run_treeshrink.py -t 5168_testPass.fasta.treefile -q "0.05 0.10" > 5168.log
```

---

I wrote a code to output the statistics from a FASTA file. Those statistics are simple and basic, which include average, median, longest, and the shortest sequence length in a FASTA file. The reason I did this is to see whether our previous set 400 bp overlapping filter is valid or not. I kinda found the relationship between the average sequence length and the quality of final tree on a fixed 400bp filter. It is obvious that the longer sequence length will be give a worse result. Instead, I modified the code to filter the minimum overlapping based on the average squence length in a FASTA file. I expected to receove better results this time.

I will test the new way, as well as incorporate the workflow I had above. Here is the command line:
```bash
while read line
do
        python remove_overlapped.py ../output/merged_ref_contigs/${line}_ref_contig_merged.fasta ${line}_reduced.fasta ${line}.nonoverlapped.fasta --num_threads 128 --filter_intensity 0.1
        mafft --preservecase --maxiterate 1000 --localpair --adjustdirection --thread 128 ${line}_reduced.fasta > ${line}_aligned.fasta
        ~/software/trimal/trimal/source/trimal -in ${line}_aligned.fasta -out ${line}_trimmed.fasta -gt 0.5
        python composition_test.py ${line}_trimmed.fasta ${line}_testPass.fasta
        python extract_filtered_seq.py ${line}_testPass.fasta ../output/merged_ref_contigs/${line}_ref_contig_merged.fasta ${line}_testPass_original.fasta
        mafft --preservecase --maxiterate 1000 --localpair --adjustdirection --thread 128 ${line}_testPass_original.fasta > ${line}_2ndaligned.fasta
        ~/software/trimal/trimal/source/trimal -in ${line}_2ndaligned.fasta -out ${line}_2ndtrimmed.fasta -gt 0.5
        python composition_test.py ${line}_2ndtrimmed.fasta ${line}_2ndtestPass.fasta
        ~/software/iqtree/iqtree-2.2.2.7-Linux/bin/iqtree2 -s ${line}_2ndtestPass.fasta -m MFP -bb 1000 -redo
done < subsample_genelist.txt
```