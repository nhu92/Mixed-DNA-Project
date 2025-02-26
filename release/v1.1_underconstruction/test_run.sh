echo "# The pipeline is recommended to run under job submission systems."
echo "# Step 0: Pipeline Clone"
git clone -n  https://github.com/nhu92/Mixed-DNA-Project.git --depth=1
cd Mixed-DNA-Project/
git checkout HEAD release/v1.1_underconstruction
git checkout HEAD release/sample_data
cd ..
mkdir MDNA_test
cd MDNA_test
cp -r ../Mixed-DNA-Project/release/v1.1_underconstruction/* ./
cp -r ../Mixed-DNA-Project/release/sample_data/* ./
rm -rf ../Mixed-DNA-Project/


# If you have your own sample to test, please ignore the following steps.
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
# ----
