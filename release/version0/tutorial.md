# Pipeline Tutorial for Plant Species Identification from Mixed Samples

Nan Hu, 2024 August

---

This tutorial provides step-by-step instructions to run a pipeline designed for identifying plant species from mixed DNA samples using the Angiosperms353 target sequencing kit and the HybPiper workflow. The pipeline is efficient and cost-effective, making it a valuable tool in various scientific and practical applications.

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
  - `HybPiper`, `fastp`, `mafft`, `fasttree`
  - Required Python libraries: `pandas`, `argparse`, `scipy`, `scikit-learn`, `numpy`, `biopython`
  `pip install numpy pandas scipy scikit-learn biopython argparse`

- **Data**:
  - Paired-end reads from mixed plant DNA samples.
  - Reference database of Angiosperms353 sequences.
  > Reference sequences name should follow Order_Family_Genus_Species format

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
python 01_exons_assembly.py -t <threads> -r1 <read1.fastq> -r2 <read2.fastq> -m <mega353.fasta> -p <project_name> -l <log_file>
```

**Arguments**:

- `-t` or `--threads`: Number of CPU threads to use.
- `-r1` or `--read1`: Path to the first read file (FASTQ format).
- `-r2` or `--read2`: Path to the second read file (FASTQ format).
- `-m` or `--mega353`: Path to the MEGA353 reference file (FASTA format).
- `-p` or `--project_name`: Project name for output files.
- `-l` or `--log_file`: Path to the log file.

### Step 2: Exon Tree Creation (`02_exon_trees.py`)

In this step, the pipeline extracts exons from the assembled sequences and builds exon trees for phylogenetic inference.

#### Usage

```bash
python 02_exon_trees.py -g <gene_list.txt> -o <overlapping_rate> -p <project_name> -l <log_file>
```

**Arguments**:

- `-g` or `--gene_list`: Path to the list of gene names.
- `-o` or `--overlapping_rate`: Rate of overlap for exon extraction.
- `-p` or `--project_name`: Project name for output files.
- `-l` or `--log_file`: Path to the log file.

### Step 3: Distance Matrix Calculation (`03_distance_matrices.py`)

This step calculates distance matrices based on the exon trees created in the previous step. These matrices are used in phylogenetic analysis to infer relationships among the taxa in the mixed sample.

#### Usage

```bash
python 03_distance_matrices.py -i <input_file> -o <output_directory> -m <method>
```

**Arguments**:

- `-i` or `--input_file`: Path to the input file containing exon trees.
- `-o` or `--output_directory`: Directory where the distance matrices will be saved.
- `-m` or `--method`: Method for distance calculation (e.g., `euclidean`, `hamming`).

### Step 4: Prediction and Identification (`04_prediction.py`)

The final step involves predicting and identifying the plant species present in the mixed sample by applying phylogenetic inference on the computed distance matrices.

#### Usage

```bash
python 04_prediction.py -i <input_distance_matrix> -o <output_file> -tl <taxonomic_level>
```

**Arguments**:

- `-i` or `--input_distance_matrix`: Path to the input distance matrix file.
- `-o` or `--output_file`: Path to save the prediction results.
- `-tl` or `--taxonomic_level`: Taxonomic level to process (choices: `o` for Order, `f` for Family, `g` for Genus, `s` for Species).

### Example Execution

To run the entire pipeline, execute each script in the following order:

```bash
# Step 1: Sequence Assembly
python 01_exon_assembly.py -t 64 -r1 01x02x03.R1.fastq -r2 01x02x03.R2.fastq -p z010203 -g gene.list.txt

# Step 2: Exon Tree Creation
python 02_exons_phylo.py -t 64 -p z010203 

# Step 3: Distance Matrix Calculation
python 03_distance_matrices.py -t 64 -p z010203 --threshold 1

# Step 4: Prediction and Identification into Order
python 04_prediction.py -i distance_matrices/matrix.csv -o predictions.csv -tl o
```

The results should predict Roaceae only for a species mixed in 3.
