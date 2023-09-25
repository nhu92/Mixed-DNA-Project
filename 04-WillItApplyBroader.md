# Will the pipeline good enough to be applied to a broader utilization?

Testing if the pipeline is suitable with different input. We expect they could still maintain at least family-level resolution. There will be some tunings but no large modifying of the current one. If it is not working, I have to develop it all over again taking the assemble genes idea.

---

As mentioned in the last chapter, I want to do two things to check the ability of our pipeline. The first is to run the pipeline to the unknown mix. Let's start with running HybPiper.

```bash
# similar command lines but do not need to run the while loop
#!/bin/bash
#SBATCH -J hyb_asbl_known_mix
#SBATCH -p nocona
#SBATCH -o log/%x.out
#SBATCH -e log/%x.err
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=3G

hybpiper assemble -t_dna ../raw/mega353.fasta \
        -r ../raw/SSBseq003.trimmed.R1.fastq.gz ../raw/SSBseq002.trimmed.R2.fastq.gz \
        --prefix unknownmix \
        --bwa \
        --cpu 32 \
        -o ../hyb_output
```

It takes 15 min to run. And then I threw it into the pipeline. The phylogeny does not have so many hits left after filtering. However, I found some patterns that suggest the taxonomy information of the unknown mix. Most of the exons seems to come from some Fabids and Malvids species, and may have some of the outgroups (Monocots, basal dicot?). I am now downloading reads  from Kew Garden database that select from the suspicious clade and them run the pipeline again to anchor these species.

The idea is to select one monocot (maybe Poaceae species?), one COM Fabids (Malpighiaceae? I was working on Salicaceae but they are all trees...), one Nitrogen fixer Fabids (Rosaceae?), one Myrtaceae, one Malvaceae, and one outgroup (Asteraceae?). In total, it is 6 species.

I will choose the latest one from each family and only use the PAFTOL source data. Also, I will keep an eye on the number of genes recovered. The material I would like to keep in consistent with silica-dried.

```
Poaceae: Gigantochloa atter ERR7621555
COM Fabids, Irvingiaceae: Irvingia gabonensis ERR4180054
N fixer Fabids, Urticaceae: Boehmeria ramiflora ERR7622000
Myrtaceae: Corymbia ficifolia ERR5034279
Malvaceae: Hildegardia barteri ERR7622285 (Herbarium)
Asteraceae: Bethencourtia palmensis ERR9230212
```

Don't forget to trim the reads!
```bash
while read line
do
fastp -i ../raw/${line}_R1.fastq.gz -I ../raw/${line}_R2.fastq.gz -o ../raw/${line}_trimmed_R1.fastq.gz -O ../raw/${line}_trimmed_R2.fastq.gz -j ../output/${line}.json -h ../output/${line}.html
done < namelist.txt
```
Submitted all the jobs. Wait and probably submit the module 2 tonight.

---

The whole exon pipeline is running for the 6 species. I created the expected topology for the comparison. 

(Gigantochloa_atter, (((Irvingia_gabonensis, Boehmeria_ramiflora), (Corymbia_ficifolia, Hildegardia_barteri)), Bethencourtia_palmensis));

Now, I am going to deal with the unknown mix data. The unknown mix seems to have low numbers of assembled contigs. I will try to lower the criteria in HybPiper assemble process.

Asteraceae: Bethencourtia palmensis ERR9230212 is not working (for unknown reason, I will shift to Monarrhenus pinifolius). Here is the updated Species list and the expected tree:
```
Poaceae: Gigantochloa atter ERR7621555
COM Fabids, Irvingiaceae: Irvingia gabonensis ERR4180054
N fixer Fabids, Urticaceae: Boehmeria ramiflora ERR7622000
Myrtaceae: Corymbia ficifolia ERR5034279
Malvaceae: Hildegardia barteri ERR7622285 (Herbarium)
Asteraceae: Monarrhenus pinifolius ERR7621651

(Gigantochloa_atter, (((Irvingia_gabonensis, Boehmeria_ramiflora), (Corymbia_ficifolia, Hildegardia_barteri)), Monarrhenus pinifolius));
```

The results seem to be overfiltered. I modified the criterion of filtering process. -i 0.8 -> 0.2, -p 0.1 -> 0.4. If this is still too stringent I will also change the critical value of the composition test from 0.05 to 0.01 to rescue some exons from merged data.

---

Back to work on Friday. D20 =18 so ... let's make some real progress to wrap up this week! I have 3 things to try on today:
1. Review the results from the pipeline.
2. Redo hybpiper for unknown mix, then merge it with 8 speices and 6 species selected tree.
3. Find a way to cut the outlier long branches from the tree. Test for the criteria.

Got a bad result from the exon pipeline of 6 species. None of the gene trees share the same topology as expected. I wonder what's going on there. I am running a full contig pipeline now.

I started the rework of hybpiper unknown mix with an additional parameter `--cov_cutoff 4`.

---

The module 2 is extremely affected by large input files (merged exons/contigs). I need to add one seletion step to ensure the robustness of the whole pipeline.

I will offer a better pipeline to get this problem solved in line automatically.

