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

Poaceae: Gigantochloa atter ERR7621555
COM Fabids, Irvingiaceae: Irvingia gabonensis ERR4180054
N fixer Fabids, Urticaceae: Boehmeria ramiflora ERR7622000
Myrtaceae: Corymbia ficifolia ERR5034279
Malvaceae: Hildegardia barteri ERR7622285 (Herbarium)
Asteraceae: Bethencourtia palmensis ERR9230212




