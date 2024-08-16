# Pipeline Tutorial for Plant Species Identification from Mixed Samples

Nan Hu, 2024 August

---

This tutorial provides step-by-step instructions to run a pipeline designed for identifying plant species from mixed DNA samples using the Angiosperms353 target sequencing kit and the HybPiper workflow. The pipeline is efficient and cost-effective, making it a valuable tool in various scientific and practical applications.

### Quick Example Run

To run the entire pipeline, execute each script in the following order:

```bash
# Step 0: Pipeline Clone
git clone -n  https://github.com/nhu92/Mixed-DNA-Project.git --depth=1
cd Mixed-DNA-Project/
git checkout HEAD release/version0
git checkout HEAD release/sample_data
cd ..
cp -r Mixed-DNA-Project/release/version0/* ./
cp -r Mixed-DNA-Project/release/sample_data/* ./
rm -rf Mixed-DNA-Project/
gunzip angiosperms353_v2_interim_targetfile.fasta.gz

  ## In-silico Mix Generation
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

# Step 1: Sequence Assembly
python 01_exon_assembly.py -t 64 -r1 01x02x03.R1.fastq -r2 01x02x03.R2.fastq -p z010203 -g gene.list.txt

# Step 2: Exon Tree Creation
python 02_exons_phylo.py -t 64 -p z010203 

# Step 3: Distance Matrix Calculation
python 03_distance_matrices.py -t 64 -p z010203 --threshold 1

# Step 4: Prediction and Identification into Order
python 04_prediction.py -i distance_matrices/matrix.csv -o predictions.csv -tl o

# Additional steps
# If want to predicted into family level, a customized reference file should be generated through:
while read line; do python pick_match_list.py ref_871/${line} ref/${line} species_sp.txt; done < gene.list.txt
# where species_sp.txt is a list of keywords in species taxonomy groups that we want to use as reference in the next round. For example, if result is Rosales, the list can just have "Rosales" in the list file.
# Then, just apply Step 2-4 to generate a prediction in family level. 
```

The results should predict **Rosales** only for DNA of species mixed in 3. If go further into family, it should predict **Rosaceae**, **Ulmaceae**, and **Canabaceae**.

## Background

Plant species identification from mixed samples is critical in fields such as ecology, conservation, and agriculture. Traditional methods face challenges due to high costs, computational inefficiencies, and incomplete databases. This pipeline addresses these issues by leveraging a novel approach that combines target sequencing and phylogenetic inference.

## 1. Data Processing

- **Input**: The process starts with raw mixed plant samples.
- **Steps**:
  - **DNA Extraction**: DNA is extracted from the mixed plant sample following established protocols (e.g., Hale et al., 2019).
  - **Target Sequencing Library Preparation**: The extracted DNA undergoes library preparation for target sequencing (as Slimp et al., 2021).
  - **Illumina Sequencing**: The prepared libraries are sequenced using Illumina technology to produce short reads, which serve as the input for the next stage.

## 2. Target Assembly

- **Input**: Short reads generated from the Illumina sequencing step.
- **Steps**:
  - **FastP Trimming & QC**: The short reads undergo quality control and trimming using FastP to produce clean reads.
  - **HybPiper**: The clean reads are assembled using HybPiper, generating BAM files, contigs in FASTA format, and exons in FASTA format. These outputs are critical for the subsequent phylogenetic inference.

## 3. Phylogenetic Inference

- **Input**: Contigs and exon FASTA files from the Target Assembly stage.
- **Steps**:
  - **MAFFT Alignment**: The contigs and exons are aligned using MAFFT against a reference panel.
  - **Trimal**: The aligned exons are trimmed to ensure consistency and remove low-quality alignments.
  - **FastTree/IQ-TREE**: The trimmed alignments are used to infer phylogenetic trees, which illustrate the evolutionary relationships among the sequences.

## 4. Prediction

- **Input**: Exon trees generated during Phylogenetic Inference.
- **Steps**:
  - **Genetic Distance Calculation**: Distance matrices are calculated from the phylogenetic trees to quantify the genetic differences among the taxa.
  - **ACS Evaluation**: An ACS (Adaptive Cluster Sampling) array is created based on the distance matrices to evaluate the species' relationships.
  - **Species Prediction**: The final species prediction is made using the ACS array, identifying the plant species present in the mixed sample.

## Prerequisites

Before running the pipeline, ensure that you have the following:

- **Software and Tools**:
  - Python 3.8+
  - `HybPiper`, `fastp`, `mafft`, `fasttree`, `seqkit`. Suggest using Conda to install
  - Required Python libraries: `pandas`, `argparse`, `scipy`, `scikit-learn`, `numpy`, `biopython`

```bash
  # Create environment for HybPiper
  conda config --add channels defaults
  conda config --add channels bioconda
  conda config --add channels conda-forge
  conda create -n pacificod hybpiper
  conda activate pacificod
  # Install dependencies
  conda install seqkit
  conda install fasttree
  conda install fastp
  pip install numpy pandas scipy scikit-learn biopython argparse
```

- **Data**:
  - Paired-end reads from mixed plant DNA samples.
  - Reference database of Angiosperms353 sequences.
  > Reference sequences name should follow Order_Family_Genus_Species format for prediction. For example: >Rosales_Rosaceae_Rose_rosa

## Pipeline Overview

The pipeline consists of four main steps:

1. **Sequence Assembly**
2. **Exon Tree Creation**
3. **Distance Matrix Calculation**
4. **Prediction and Identification**

### Step 1: Sequence Assembly (`01_exons_assembly.py`)

This step involves trimming and assembling short paired-end reads using `fastp` and `HybPiper`.

#### Usage

```bash
python 01_exons_assembly.py -t <threads> -r1 <read1.fastq> -r2 <read2.fastq> -m <mega353.fasta> -p <project_name> -g <gene_list>
```

**Arguments**:

- `-t` or `--threads`: Number of CPU threads to use.
- `-r1` or `--read1`: Path to the first read file (FASTQ format).
- `-r2` or `--read2`: Path to the second read file (FASTQ format).
- `-m` or `--mega353`: Path to the MEGA353 reference file (FASTA format). Default is `angiosperms353_v2_interim_targetfile.fasta`.
- `-p` or `--project_name`: Project name for output files.
- `-g` or `--gene_list`: Path to the selected gene list.
- `-ov` or `--overlapping_rate`: Overlapping to consider the same exon (0-1, default 0.8).

### Step 2: Exon Tree Creation (`02_exon_trees.py`)

In this step, the pipeline extracts exons from the assembled sequences and builds exon trees for phylogenetic inference.

#### Usage

```bash
python 02_exon_trees.py -t <threads> -e <exon_dir> -r <ref_dir> -p <project_name> -g <gene_list>
```

**Arguments**:

- `-t` or `--threads`: Number of CPU threads to use.
- `-e` or `--input_exon`: Directory of extracted exon sequences. Default is `02_exon_extracted`.
- `-r` or `--input_exon`: Directory of reference alignments. Default is `ref`.
- `-p` or `--project_name`: Project name for output files.
- `-g` or `--gene_list`: Path to the list of gene names. Default is generated from Step 1 `gene_list.txt`.

### Step 3: Distance Matrix Calculation (`03_distance_matrices.py`)

This step calculates distance matrices based on the exon trees created in the previous step. These matrices are used in phylogenetic analysis to infer relationships among the taxa in the mixed sample.

#### Usage

```bash
python 03_distance_matrices.py -t <threads> -p <project_name> -g <gene_list> --threshold <float_num> [--use_flag] 
```

**Arguments**:

- `-t` or `--threads`: Number of CPU threads to use.
- `-p` or `--project_name`: Project name for output files.
- `-g` or `--gene_list`: Path to the list of gene names. Default is generated from Step 1 `gene_list.txt`.
- `--threshold`: Parameter to control specificity of candidate branches on phylogeny. Default is `1.96`, which means only genetic similarity is over `1.96*SD+Mean` will be considered.
- `--use_flag`: A boolean parameter to mark if all the non-outliers distances set to 0.

### Step 4: Prediction and Identification (`04_prediction.py`)

The final step involves predicting and identifying the plant species present in the mixed sample by applying phylogenetic inference on the computed distance matrices.

#### Usage

```bash
python 04_prediction.py -i <input_ACS_file> -o <output_file> -tl <taxonomic_level>
```

**Arguments**:

- `-i` or `--input_distance_matrix`: Path to the input cumulative distances file (*cumulative_dist.csv).
- `-o` or `--output_file`: Path to save the prediction results.
- `-tl` or `--taxonomic_level`: Taxonomic level to process (choices: `o` for Order, `f` for Family, `g` for Genus, `s` for Species).
