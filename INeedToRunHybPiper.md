# Wish I could rerun the HybPiper!
This is a diary-style markdown for documenting steps for a rerun of HybPiper on Mixed DNA samples

---

1. Delete conda, re-install conda

I had to delete conda in the HPCC of TTU since it became super slow for the past couple of years. Now I am graduated and temporarily free from worrying about reproducing my old analysis (I hope) (And I am working on documenting those pipelines, slowly. See [HomologousAnnotationPipeline](https://github.com/gudusanjiao/HomologousAnnotationPipeline)). As suggested on Anaconda website, there are two ways to uninstall conda from remote computer clusters - the full uninstall and the 'nearly' full uninstall. The full uninstall is automatically operated by conda but it requires a functional conda (which I don't have and that's why I want to reinstall right?). Thus, I did a manual uninstall, which is basically just deleting everything from conda directories and wiping out most of the traces from user files.
As suggested, here are command lines and additional steps to take into action.
```bash
rm -rf ~/anaconda3
rm -rf ~/conda
```
Also, search for any "suspicious" files created by conda in `/home/<user>`.

After this, search for files that look like `.bash*` (for example: `.bashrc`). It requires to use `ls -a` to show those hidden files. Delete anything that popped up in those files that you think relates to conda.

Now, I get a fresh and clean "conda-free" environment for my HPCC account. It's about time to restart everything!

OK! I can't wait to reinstall conda now! To be safe (Avoid issues I faced before) and fast, I reinstalled the lite version of conda - [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/index.html). Choose the correct one (normally the first one of the Linux platform works), then download it.
```bash
wget TheAddressOfThePackage.sh
sh TheFileNameOfTheThingYouJustDownloaded.sh
```
Then, hit 'q', 'yes', and 'yes' to finish the installation. Restart the connection to HPCC and type `conda --help` to check the availability of conda.

Do something to prevent the future collapse of conda! Replace the original solver of conda with mamba
```bash
# why this step also takes forever to finish????
# Update: It takes 45min...
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
```
Starting today, organize the conda environment and separately install different packages in different environments.
```bash
conda create -n hybpiper
# visit new env by:
conda activate hybpiper
# leave the env by:
conda deactivate
```

Install hybpiper via conda:
```bash
conda install hybpiper
# failed with missing libgcc-ng >=12, I am about to install it with original solver.
# explore issue report online
# finally solved by using:
conda create -n hybpiper -c chrisjackson-pellicle -c bioconda -c conda-forge hybpiper
```
Test if it is installed correctly!
```bash
conda activate hybpiper
hybpiper --help
hybpiper check_dependencies
```

2. Prepare the input file for HybPiper

The HybPiper requires two input files: a target sequence FASTA (which is mega353.fasta from Angiosperm353) and a sequence capture reads file. The reads need to be trimmed and passed through quality control using Trimmomatic and fastp. I directly used data from Courtney so I can skip this step (need to check if errors appear in the future).

According to the notes from Courtney, I assembled all trimmed reads using `hybpiper assemble` with parallel attempts.

```bash
#!/bin/bash
#SBATCH -J hyb_asbl_REPLACE
#SBATCH -p nocona
#SBATCH -o log/%x.out
#SBATCH -e log/%x.err
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=3G

hybpiper assemble -t_dna ../raw/mega353.fasta \
        -r ../raw/REPLACE.R1.trimmed.fastq ../raw/REPLACE.R2.trimmed.fastq \
        --prefix REPLACE \
        --bwa \
        --cpu 64 \
        -o ../hyb_output
```
I need to generate a list file containing all of the species names to replace the "REPLACE" in my code.
```bash
while read line; do sed s/REPLACE/$line/g 01-hyb_assemble_REPLACE.sh | sbatch; done < namelist.txt
```

It's 5:30pm now. I leave all the scripts on screen or submitted to slurm. It's time for wrapping up today.

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

