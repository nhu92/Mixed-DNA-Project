# Using the contig contribution to estimate the OTUs and other function to adjust the ACS

I will record the information that each exon was came from which contig to calculate the contribution of the contigs in each taxonomy group. Also, this developing note will track the optimizing progress of the pipeline to output better predictions.

---

To estimate the OTUs (Operational Taxonomic Units), there are several thoughts in my mind. QIIME uses similarity level at 97% to consider the same OTUs, however, not quite applicable in our data. The reason is that the plant references are still very sparse to really narrow down to the genus/species. Also, our solution of using the hierarchical approach would have different similarity within various taxonomy groups. Thus, I would like to trust HybPiper's ability to assemble different origins of DNA from the mixed samples. Whether the same taxonomy group has majority multiple contigs from the same gene to contribute into, will become the indicator that the input data has multiple species within the same group. I will go through this logic to test the ability of our pipeline.

The downside of this approach can be: 1) This assembled contigs are shorter than the full length exons. Thus, there will be more than one copy for each exons on a gene. They follow a distribution, while, not really equal to the number of species in the mix. We need an algorithm to estimate the real number of species from this information; 2) Plants are naturally favor of polyploidy. Our pipeline is very hard to distinguish the difference between multiple species in the same taxonomy group or it is polyploidy (especially for allopolyploidy). I do not have any potential polyploidy samples to test this but I guess there will be differences in the contig/exon distribution between two independent species in the same group and the same species with two copies of genes.

---

I fixed a critical issue from the pipeline. The previous `run_command()` function will quit the entire pipeline when there are any errors happened on the command it takes. Thus, when the `ls` is failed to take the list of exons from the gene (mainly because HybPiper does not generate a useful exonerate output), the entire pipeline will crash and exit. The solution is to take a `run_command` function into "critical" or "non-critical" with a parameter to handle. For the situation like the `ls` command, it is not critical but for some steps (e.g. HybPiper), it should be critical. I will change all the code to avoid this issue and rerun the testing pipeline on the family level.
