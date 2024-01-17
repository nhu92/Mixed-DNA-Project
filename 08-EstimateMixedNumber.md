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