# Wish I could rerun the HybPiper!
This is a note-style markdown for documenting steps for a rerun of HybPiper on Mixed DNA samples

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
#SBATCH -J hyb_asbl
#SBATCH -p nocona
#SBATCH -o log/hyb_asbl%x.out
#SBATCH -e log/hyb_asbl%x.err
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=3G

hybpiper assemble -t_dna ../raw/mega353.fasta \
        -r REPLACE.R*.trimmed.fasta \
        --prefix REPLACE \
        --bwa \
        --cpu 64 \
        --run_intronerate \
        --start_from exonerate_contigs
```
I need to generate a list file containing all of the species names to replace the "REPLACE" in my code.
```bash
while read line; do sed 's/REPLACE/$line/g' 01-hyb_assemble_REPLACE.sh | sbatch; done < namelist.txt
```

It's 5:30pm now. I leave all the scripts on screen or submitted to slurm. It's time for wrapping up today.

