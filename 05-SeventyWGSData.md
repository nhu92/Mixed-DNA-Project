# Using well developed genome sequencing to build a reference gene list
When we have the gene sequence in a good quality, it will be very easy to use mafft --add to align any unknown sequences within a reasonable computing time.

---

Let me generate a species tree of these 70 species? Manually? I need a taxonomy table for them first. 

I checked gene 4471 (the first one). It seems that all of the exon extractions have some hits with the established alignment from the 70 species 4471. However, for long sequences, not the entire exon extracts can be mapped to the alignment. Some fragments performs better than others. The tree is fine right now with tested species (Salvia) grouped with monkey flower (Mimulus ringens) for most of times. This may indicate a distance method would be possible to identify the species.

My current plan is, from exonerate stats file, extract each hit and make phylogeny independently to assess the phylogeny positions. For each tree, calculate the distance of the targeted species to other species and list them out.

Here is a basic thought about this code:
```python
# extract ranges from exonerate stat file, column 6
    ## recognize the structure [(a,b),(c,d),(e,f)]
# check overlaps, then name each range
    ## 90% exon range overlap will be considered as the same exon
# pull the exons from contig file and named by the exon info, column 13
    ## recognize the structure [(a,b),(c,d),(e,f)]
# output the exon with the exon orders, this output will be the input file for mafft --add
```

