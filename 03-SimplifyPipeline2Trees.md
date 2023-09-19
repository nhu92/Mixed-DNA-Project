# An easier way to filter the data. Use HybPiper again!

Utilize the exonerate data to extract conservative region across all the hits contigs.
---

As the original pipeline works, the remaining issue now is to expand the "in" pool of gene tree construction. Rather than adjust the filtering criteria, the more valid way is to direct extract the targeted gene regions from the output tsv from the exonerate run of HybPiper.

Let's check the format of the exonerate output. The file is stored in `{hybpiper_dir}/{gene_name}/{hybpiper_proj_name}/exonerate_stats.tsv`.

Here is a preview of this file:

|                                       | query_id  | query_length | hit_id                            | query_HSP_range_limits_original | query_HSP_range_limits_trimmed | query_HSPFragment_ranges                                                       | hit_percent_similarity_original | hit_percent_similarity_trimmed | hit_strand | hit_HSP_range_limits_original | hit_HSP_range_limits_trimmed | hit_HSPFragment_ranges_original                                                              | hit_HSPFragment_ranges_trimmed                                                               | 3-prime_bases_trimmed |
|---------------------------------------|-----------|--------------|-----------------------------------|---------------------------------|--------------------------------|--------------------------------------------------------------------------------|---------------------------------|--------------------------------|------------|-------------------------------|------------------------------|----------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|-----------------------|
| Hits filtered > 55 percent similarity | VLNB-5974 | 234          | NODE_8_length_1347_cov_20.507377  | (0, 132)                        | (0, 132)                       | [(0, 5), (5, 21), (21, 34), (35, 68), (68, 132)]                               | 92.48                           | 92.48                          | 1          | (36, 1303)                    | (36, 1303)                   | [(36, 51), (290, 338), (468, 508), (603, 704), (1111, 1303)]                                 | [(36, 51), (290, 338), (468, 508), (603, 704), (1111, 1303)]                                 | N/A                   |
|                                       | VLNB-5974 | 234          | NODE_7_length_1600_cov_29.737950  | (2, 234)                        | (2, 234)                       | [(2, 21), (21, 34), (35, 68), (68, 132), (133, 165), (166, 204), (205,   234)] | 74.15                           | 74.15                          | -1         | (199, 1530)                   | (199, 1530)                  | [(1473, 1530), (1313, 1353), (1127, 1228), (770, 963), (594, 693), (395,   513), (199, 284)] | [(1473, 1530), (1313, 1353), (1127, 1228), (770, 963), (594, 693), (395,   513), (199, 284)] | N/A                   |
|                                       | VLNB-5974 | 234          | NODE_10_length_1012_cov_73.506215 | (0, 68)                         | (4, 68)                        | [(0, 21), (21, 34), (35, 68)]                                                  | 81.16                           | 86.15                          | 1          | (159, 546)                    | (171, 546)                   | [(159, 222), (303, 343), (445, 546)]                                                         | [(171, 222), (303, 343), (445, 546)]                                                         | N/A                   |

The useful information is the name of the sequence (Column 4) and the hits match range (Column 12). We will use these information to extract the contig that has the hit and trim the contig to the hit only segment.

Have a code to extract the exon from the fasta and output. The code is tested and works fine:

```bash
python exon_hits_collect.py ../hyb_output/knownmix/ ../output/exon_hits/
```

There are some blank sequences in the FASTA. What happened? (Deal with it tomorrow, really affect the automatic process).
> Update: 09/19/2023 Fixed. The issue was caused by the range of selection. The export of exonerate gives a range can start from 0 but in our script we select [start - 1: end] where when we hit a range like (0, 599) it will extract nothing. A simple fix is taken (and you should know how to :))

The filtered one looks nice. I neede to generate a mixed tree out of it and see the resolution.

---

I want to separate different exons as different genes, then extract the segment that mapped to the exons.

I would like to edit the remove_nonoverlapped.py to allow a parameter about the size selection intensity. Also, I want to rename the proportion of overlapped parameter.

