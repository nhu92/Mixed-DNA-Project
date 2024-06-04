# Let's make some nice plots for the poster.

I will unify the plot generation using R::ggplot2. 

---
Here is the note from the meeting with Dr. Johnson.

Notes for Poster
1. Use of 800 species as flowering plant phylogeny reference
2. Development of Z-score method for assessing similarity across many loci
    a. topo-map figure, example tree(s)
    b. classifier logic (using topology and the MRCA support values)
3. Positive/Negative predictive value analysis (include definitions)
    a. how many reference samples we need to ID to order: number of genes and number/phylogeneitc spread of taxa
    b.mixture numbers versus accuracy - less “busy” figure by just choosing mixes of 3, 6, and 11
4. Specificity/Precision analysis (minimum reads to ID)
5. Runtime analysis

The first figure would be the regeneration of the topo map. This figure was removed from the current analyzing pipeline. I will recall the code and try to add some reasonable edition to it.

This figure requires the astral tree of the testing reference. Since we use the new A353v2, I would like to evaluate the tree before the plot generation. I could reuse the pipeline from 02-reference.sh without generate the contour map. The pipeline should finish before the exon tree generation process.

Find some old code file and pipelines to do this. I may remove the SD for a reason.

---

I will also generate a contour map for the result. I had never used that plot for the 100 species map. I think it might be messy but I will check.

Glad to pick up some old scripts generated before. I almost forget them already. The contour map will be manually generated from the two inputs: the reference PCoA result and the cumulative similarity result. I might need to fix the PCoA result by making the seq name into specific hierarchies. For example, we need the level at orders/families. Also, the final contour map does not like a very long descriptional naming.

This is the small pipeline to generate the expected contour map. We do not need to run the ref_panel pipeline every time for this work. If the reference file does not change (which means, both the number of genes and the reference species), we can use the same `_refPCoA.csv` as the input file.

```bash
proj_name=
level=
python categ_PCoA.py 113_30_refPCoA.csv 113_30_refPCoA_${level}.csv ${level}
dos2unix 113_30_refPCoA_${level}.csv
dos2unix ${proj_name}.Zaxis.csv
paste -d, 113_30_refPCoA_${level}.csv ${proj_name}.Zaxis.csv > ${proj_name}_3dcoord${level}.csv
python contour_optima_normalized.py s${proj_name}_3dcoord${level}.csv ${proj_name}_contourmap${level}.svg 0
```

---

I will start to rerun the 100 tests from the Kew Garden database. Should I remove the individiuals that we found some issues with? Or I can replace those species with other ones from the Kew? Also, I wish to sample each species with even times to avoid bias. Also, I want to generate a super mix with reduced reads number.

I will check each species reads number, then, I will thin down the read number for each species into 150k reads to avoid a huge input file. 

---

I need to generate some example trees for the display. It is better to color them in a good way. I probably need to show the heatmap with the topology on side. 

---

This is the process of generating new sets of 100 mixes.

I will generate 40 3x mixes (120sp, 6/sp), 40 6x mixes (240sp, 12/sp), and 20 10x mixes (200sp, 10/sp). The detailed processes is divided below:

1. Thin down the reads into 1000k (very sufficient) for each input.
2. Tablize input name to make mix, each column is a random redistribution of 20 names
3. Merge read according to the random table.

The code for making the mix. 
```R
namelist <- c("01x", "02x", "03x", "04x", "07x", "08x", "09x", "10x", "12x", "13x", "15x", "16x", "18x", "19x", "22x", "24x", "25x", "26x", "28x", "30x")
mix3 <- t(replicate(40, sample(namelist, 3, replace = FALSE)))
mix6 <- t(replicate(20, sample(namelist, 6, replace = FALSE)))
mix10 <- t(replicate(20, sample(namelist, 10, replace = FALSE)))

```

---

Submitted 100 jobs on 10 genes tests. I expect some of them failed due to limited gene number. I will expand those tests to more genes tomorrow.

---

This is the new evaluation pipeline for the new 100 tests:
```bash
while read line; do python cumu_order.py ../gene10/${line}/${line}.cumulative_dist.csv ../eval_gene10/${line}.cumulative_dist_order.csv; python evaluate.py ../eval_gene10/${line}.cumulative_dist_order.csv candidates_order.tsv ../eval_gene10/${line}.order_eval.csv; done < mix_namelist.txt
```

---

Back to work at TTU today! I gonna check the parallel script to make sure the parallel is perfectly working on the distance matrix.

The previous parallel was running with 4 jobs, and the process was executed by each gene. The speed was limited both by the number of jobs and the number of exons in each gene. The whole process saved around 1 hour (from 5 to 4). I will try increase the job number in the parallel to see if it drastically reduces the run time. Also, I will work on a script to "flatten" the input tree without categorize them within genes. If we could process all the exon trees at the same time, I believe it will save more time.

The reason I do so for the parallel is that we found the accuracy will be obviously increased by using all the references from the A353v2.

The 64 threads run finished the process less than 30 min, which is a very strong improvement from the original pipeline. I will use this to rerun the 800+ references test.

The current version I am using:
```bash
#!/bin/bash
#SBATCH -J parallel_test
#SBATCH -p nocona
#SBATCH -o log/%x.out
#SBATCH -e log/%x.err
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=3G

# This is a pipeline for analyzing raw reads from a mixed sample to the identification.
# This pipeline is to treat with the mixed reads. There is another pipeline to generate the reference panel.

# Declare the constants
threads=64
read1=01x02x13.R1.fastq.gz
read2=01x02x13.R2.fastq.gz
mega353=mega353.fasta
proj_name=01x02x13
ref_alignment=ref # The directory of reference alignment

tree_dir=all_trees_par

# Read each gene name from gene_list.txt and process files for each gene in parallel
cat ./gene_list.txt | parallel --jobs 64 '
    gene_name_shorter={};
    ls "./all_trees_par/${gene_name_shorter}"*"tre" > ./loop.treelist.txt;
    i=1;
    while read filename; do
        python matrix_ult.py -t "${filename}" -n "./all_trees/${gene_name_shorter}.${i}.list.txt" -o "./all_trees_par/${gene_name_shorter}.${i}.matrix";
        cp "./all_trees_par/${gene_name_shorter}.${i}.matrix" "./all_trees_par/${gene_name_shorter}.${i}.cleaned.csv";
        ((i++));
    done < ./loop.treelist.txt;
'

rm ./loop.treelist.txt

```


The better version I gonna test:

```bash
# flatten the read in treefiles
tree_dir=all_trees_par

for file in ${tree_dir}/*_exon_*.tre; do
  # Extract the parts of the filename
  base=${file%%_exon_*}
  exon=${file#*_exon_}
  exon=${exon%%.*}
  # Construct the new filename
  newname="${base}.${exon}.tre"
  # Rename the file
  mv "${tree_dir}/$file" "${tree_dir}/$newname"
done

# the parallel part
###### ADD a statement to create a file list (in progress)
ls ${tree_dir}/*.tre | sed 's/.tre//g' > exon_list.txt
cat ./exon_list.txt | parallel --jobs 64 '
    filename={};
    python matrix_ult.py -t "${filename}.tre" -n ./all_trees/${filename}.list.txt -o "./all_trees_par/${filename}.matrix";
    cp "./all_trees_par/${filename}.matrix" "./all_trees_par/${filename}.cleaned.csv";
'

```

---

I am working on the poster for IBC XX. I'll gradually update te progress. The first thing is to submit a job array for 150 gene test. I already optimized the pipeline for a parallel run. Thus, the runtime would be reduced a lot. 

---

Finished most of the plotting code. The 150 gene test runs not as expected. The pipeline generated a lot of empty files during exon extractions. I tried to figure out the issue. Does it come with the blank line of the `gene.list.txt`? I will double check if the pipeline works properly.

Figured out! It is the blank line issue! The blank line will be read into a `find` command for `*${gene_name}*` in the pipeline. The empty line will replace the `${}` to produce `**` in finding line, which will generate numerous files. And the whole process was produced twice in a multiplied way.

I am now rethinking about the confusion matrix generated from out pipeline. Currently, the confusion matrix is generated from counting the 4 parameters in each run, which is not that correct. There is a more "standard" way to do this in my mind. For each species in the pool, split tests into two groups: tests with the target species and the tests without the species. Then, we will evaluate if the species is predicted in each groups as TP and FP, as well as the non-predicted counts as FN and TN. The detailed table is below. Finally, summarize all the species to make a averaged predicting value across all the threshold choice.

|                         | With Species A Entries | Species A Free Entries |
|-------------------------|------------------------|------------------------|
| Species A Predicted     | True Positives         | False Positives        |
| Species A Not Predicted | False Negatives        | True Negatives         |

---

I just realized that, only using the prepared species to estimate the data is not sufficient. All of the candidate orders should be evaluated, even most of them do not show up in the sampling. They could contribute to the FP and TN counts and contribute to the positive predictions. I will manually calculate that in this case if the parameter estimation is different from previous method.

I find the two different methods may have the same number in those parameters.

---

I reviewed the results from the gene150 and 871species run and compared them with the gene30 runs. I found the boundary of the cumulative similarity is some how blurry because the source of closed related species. We used a lot of logics already to deal with this issue but it is still affecting the results. I come up with the de-correlation methods based on a matrix treating system. I make some samples using ChatGPT. I will apply this to our pipeline tomorrow to see if anything could act better than before.

Certainly! Here is the logic outlined in markdown style:

### Logic to Remove Correlation in a Matrix

#### Step-by-Step Process

1. **Define the Original Matrices**
   - The original correlation matrix and the matrix with \( Q \).

2. **Conversion of Matrices**
   - Convert each element of the initial matrices using the formula \( \frac{1}{1+x} \).

3. **Adjusting for Independence**
   - Use the converted matrices to adjust the correlation matrix such that \( Q \) becomes independent of \( a \), \( b \), and \( c \).

#### Detailed Steps

1. **Original Correlation Matrices:**
   ```plaintext
   Correlation Matrix:
   a b c
   a 0 0.1 0.9
   b 0.1 0 1
   c 0.9 1 0

   Correlation Matrix with Q:
   a b c Q
   a 0 0.1 0.9 0.4
   b 0.1 0 1 0.3
   c 0.9 1 0 0.8
   Q 0.4 0.3 0.8 0
   ```

2. **Conversion Function:**
   ```python
   def convert_matrix(matrix):
       return 1 / (1 + matrix)
   ```

3. **Conversion of Matrices:**
   ```python
   converted_new_cor_matrix = convert_matrix(new_cor_matrix)
   converted_new_cor_matrix_with_q = convert_matrix(new_cor_matrix_with_q)
   ```

4. **Function to Adjust \( Q \) to Make it Independent:**
   ```python
   def adjust_q(cor_matrix, cor_matrix_with_q):
       q_index = cor_matrix_with_q.shape[0] - 1
       
       # Placeholder for the adjusted matrix
       adjusted_matrix = cor_matrix_with_q.copy()
       
       # Iterate over the rows and columns of Q
       for i in range(q_index):
           for j in range(q_index):
               if i != j:
                   adjusted_matrix[q_index, i] -= cor_matrix[i, j] * cor_matrix_with_q[j, q_index]
                   adjusted_matrix[i, q_index] = adjusted_matrix[q_index, i]
       
       # Set the diagonal element of Q to 0 (no self-correlation)
       adjusted_matrix[q_index, q_index] = 0
       
       return adjusted_matrix
   ```

5. **Adjust the New Matrix:**
   ```python
   adjusted_new_matrix = adjust_q(converted_new_cor_matrix, converted_new_cor_matrix_with_q)
   ```

6. **Convert to DataFrame for Better Visualization:**
   ```python
   adjusted_new_df = pd.DataFrame(adjusted_new_matrix, columns=['a', 'b', 'c', 'Q'], index=['a', 'b', 'c', 'Q'])
   ```

7. **Output the Adjusted Correlation Matrix:**
   ```plaintext
   Adjusted Correlation Matrix with Q (New Data):
   a         b         c         Q
   a  1.000000  0.909091  0.526316 -0.277413
   b  0.909091  1.000000  0.500000 -0.157898
   c  0.526316  0.500000  1.000000 -0.205000
   Q -0.277413 -0.157898 -0.205000  0.000000
   ```

This logic outlines the steps taken to transform and adjust the correlation matrix, ensuring that \( Q \) is independent of \( a \), \( b \), and \( c \).