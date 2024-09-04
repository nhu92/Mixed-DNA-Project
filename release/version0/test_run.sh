echo "# The pipeline is recommended to run under job submission systems."
echo "# Step 0: Pipeline Clone"
git clone -n  https://github.com/nhu92/Mixed-DNA-Project.git --depth=1
cd Mixed-DNA-Project/
git checkout HEAD release/version0
git checkout HEAD release/sample_data
cd ..
mkdir MDNA_test
cd MDNA_test
cp -r ../Mixed-DNA-Project/release/version0/* ./
cp -r ../Mixed-DNA-Project/release/sample_data/* ./
rm -rf ../Mixed-DNA-Project/
gunzip angiosperms353_v2_interim_targetfile.fasta.gz

echo "## In-silico Mix Generation"
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR762/001/ERR7621631/ERR7621631_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR762/001/ERR7621631/ERR7621631_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR762/002/ERR7621392/ERR7621392_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR762/002/ERR7621392/ERR7621392_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR762/007/ERR7621767/ERR7621767_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR762/007/ERR7621767/ERR7621767_2.fastq.gz

seqkit sample -n 10000000 -s 100 ERR7621631_1.fastq.gz > 01_1.fastq
seqkit sample -n 10000000 -s 100 ERR7621631_2.fastq.gz > 01_2.fastq
seqkit sample -n 10000000 -s 100 ERR7621392_1.fastq.gz > 02_1.fastq
seqkit sample -n 10000000 -s 100 ERR7621392_2.fastq.gz > 02_2.fastq
seqkit sample -n 10000000 -s 100 ERR7621767_1.fastq.gz > 03_1.fastq
seqkit sample -n 10000000 -s 100 ERR7621767_2.fastq.gz > 03_2.fastq

cat 01_1.fastq 02_1.fastq 03_1.fastq > 01x02x03.R1.fastq 
cat 01_2.fastq 02_2.fastq 03_2.fastq > 01x02x03.R2.fastq 

echo "# Step 1: Sequence Assembly"
python 01_exons_assembly.py -t 64 -r1 01x02x03.R1.fastq -r2 01x02x03.R2.fastq -p z010203 -g gene.list.txt

echo "# Step 2: Exon Tree Creation"
python 02_exon_trees.py -t 64 -p z010203 

echo "# Step 3: Distance Matrix Calculation"
python 03_distance_matrices.py -t 64 -p z010203 --threshold 1

echo "# Step 4: Prediction and Identification into Order"
python 04_prediction.py -i z010203.cumulative_dist.csv -o predictions.csv -tl o

# ----------------
# Additional steps
# If want to predicted into family level, a customized reference file should be generated through:
while read line; do python pick_match_list.py ref_871/${line} ref_family/${line} species_sp.txt; done < gene.list.txt
# where species_sp.txt is a list of keywords in species taxonomy groups that we want to use as reference in the next round. For example, if result is Rosales, the list can just have "Rosales" in the list file.
# Then, just apply Step 2-4 to generate a prediction in family level. 

# Remove folders start with 03_, 04_, restart from alignment
python 02_exon_trees.py -t 64 -p z010203_fam -r ref_family
python 03_distance_matrices.py -t 64 -p z010203_fam --threshold 1
python 04_prediction.py -i z010203_fam.cumulative_dist.csv -o predictions.csv -tl f