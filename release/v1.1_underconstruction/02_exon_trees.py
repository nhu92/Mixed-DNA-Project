#!/usr/bin/env python3
"""
02_exon_trees.py – Align exon sequences and build gene trees.
This script processes exon sequences for each gene, aligning them to a reference
alignment, trimming the alignment, and constructing a phylogenetic tree.

It supports parallel processing of multiple genes and logs the status of each step.
It requires the following tools:
- MAFFT for sequence alignment
- trimAl for trimming alignments
- FastTree for phylogenetic tree construction
It can be run with command-line arguments or a configuration file.

"""
import os
import glob
import argparse
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from pipeline_utils import log_status, run_command, load_config, is_valid_project_name

def all_sequences_meet_minimum_length(fasta_path, min_length=80):
    """Check if all sequences in the FASTA file are at least `min_length` bases long."""
    for record in SeqIO.parse(fasta_path, "fasta"):
        if len(record.seq) < min_length:
            return False
    return True

def process_gene_exon_alignment(gene_name, threads, input_dir, ref_dir, output_dir, log_file, min_size):
    """
    Align and build tree for all exon contigs of a given gene.
    - Aligns contigs to the reference alignment using MAFFT.
    - Trims the alignment with trimAl.
    - Constructs a phylogenetic tree with FastTree.
    Logs progress for each exon and skips those not meeting length criteria.
    """
    # Find all exon FASTA files for this gene
    try:
        pattern = os.path.join(input_dir, f"*{gene_name}*.fasta")
        exon_files = glob.glob(pattern)
    except Exception as e:
        log_status(log_file, f"List exons for {gene_name}: FAILURE")
        print(f"Error listing exons for {gene_name}: {e}")
        return
    log_status(log_file, f"List exons for {gene_name}: SUCCESS")
    if not exon_files:
        log_status(log_file, f"No exon files found for {gene_name}, skipping.")
        return
    exon_files.sort()
    for i, exon_path in enumerate(exon_files, start=1):
        # Enforce minimum exon length
        if not all_sequences_meet_minimum_length(exon_path, min_size):
            log_status(log_file, f"Skipping {os.path.basename(exon_path)} (sequences < {min_size} bp)")
            continue
        # Alignment with MAFFT
        ref_alignment = os.path.join(ref_dir, f"{gene_name}.fasta")
        aligned_out = os.path.join(output_dir, f"{gene_name}_exon_{i}_aligned.fasta")
        mafft_cmd = (
            f"mafft --preservecase --maxiterate 1000 --localpair --adjustdirection "
            f"--thread {threads} --addfragments {exon_path} {ref_alignment} > {aligned_out}"
        )
        run_command(mafft_cmd, f"MAFFT alignment for {gene_name} exon {i}", log_file)
        # Trim alignment with trimAl
        trimmed_out = os.path.join(output_dir, f"{gene_name}_exon_{i}_trimmed.fasta")
        trimal_cmd = f"trimal -in {aligned_out} -out {trimmed_out} -gt 0.5"
        run_command(trimal_cmd, f"Trim alignment for {gene_name} exon {i}", log_file)
        # Build tree with FastTree
        tree_out = os.path.join(output_dir, f"{gene_name}_exon_{i}.tre")
        fasttree_cmd = f"fasttree -gtr -gamma -nt {trimmed_out} > {tree_out}"
        run_command(fasttree_cmd, f"Tree construction for {gene_name} exon {i}", log_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align exons and build exon trees for each gene.")
    parser.add_argument("-c", "--config", help="Path to config file (YAML/JSON/TOML)")
    parser.add_argument("-t", "--threads", type=int, help="Number of CPU threads for parallel processing")
    parser.add_argument("-e", "--input_exon", help="Directory of extracted exon FASTA files", default="02_exon_extracted")
    parser.add_argument("-r", "--ref_alignment", help="Directory of reference alignments", default="ref")
    parser.add_argument("-g", "--gene_list", help="Path to gene list file", default="gene_list.txt")
    parser.add_argument("-p", "--proj_name", help="Project name identifier")
    parser.add_argument("-o", "--output_dir", help="Directory for output trees", default="03_phylo_results")
    parser.add_argument("-m", "--min_exon_size", type=int, help="Minimum exon length to include", default=80)
    args = parser.parse_args()

    # Load config if provided
    config = {}
    if args.config:
        config = load_config(args.config)
    # Determine parameters (CLI overrides config)
    threads = args.threads if args.threads is not None else config.get('threads')
    proj_name = args.proj_name or config.get('proj_name')
    input_exon_dir = args.input_exon if args.input_exon != parser.get_default('input_exon') else config.get('input_exon', "02_exon_extracted")
    ref_dir = args.ref_alignment if args.ref_alignment != parser.get_default('ref_alignment') else config.get('ref_alignment', "ref")
    gene_list_path = args.gene_list if args.gene_list != parser.get_default('gene_list') else config.get('gene_list', "gene_list.txt")
    output_dir = args.output_dir if args.output_dir != parser.get_default('output_dir') else config.get('output_dir', "03_phylo_results")
    min_size = args.min_exon_size if args.min_exon_size != parser.get_default('min_exon_size') else config.get('min_exon_size', 80)
    if threads is None or not proj_name:
        parser.error("Required parameters missing: threads and proj_name must be specified.")
    if not is_valid_project_name(proj_name):
        parser.error(f"Project name '{proj_name}' contains invalid characters.")
    threads = int(threads)
    min_size = int(min_size)
    # Initialize log file
    log_file = f"{proj_name}_02_exons_phylo.log"
    if os.path.exists(log_file):
        os.remove(log_file)
    log_status(log_file, "Pipeline started with the following parameters:")
    log_status(log_file, f"  Threads: {threads}")
    log_status(log_file, f"  Input Exon Directory: {input_exon_dir}")
    log_status(log_file, f"  Reference Alignment Directory: {ref_dir}")
    log_status(log_file, f"  Gene List: {gene_list_path}")
    log_status(log_file, f"  Project Name: {proj_name}")
    log_status(log_file, f"  Output Directory: {output_dir}")
    log_status(log_file, f"  Minimum Exon Size: {min_size}")
    os.makedirs(output_dir, exist_ok=True)
    log_status(log_file, f"Created directory {output_dir}")
    # Read gene list and execute alignments/trees in parallel
    with open(gene_list_path, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(process_gene_exon_alignment, gene, threads, input_exon_dir, ref_dir, output_dir, log_file, min_size)
            for gene in genes
        ]
        # Ensure all threads complete (and raise exceptions if any)
        for future in futures:
            future.result()
    log_status(log_file, "Pipeline completed successfully.")
    print(f"Pipeline completed. Check {log_file} for details.")
