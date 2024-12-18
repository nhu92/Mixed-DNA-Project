import subprocess
import os
import argparse
import pandas as pd
from Bio import SeqIO
from datetime import datetime
import re
import sys

def is_valid_project_name(project_name):
    # Check if it starts with a letter, doesn't contain "NODE", and only contains letters, numbers, and underscores
    if re.match(r'^[A-Za-z0-9_-]*$', project_name) and "NODE" not in project_name:
        return True
    else:
        return False

# Function to log status with timestamp
def log_status(log_file, message):
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(log_file, 'a') as log:
        log.write(f"[{timestamp}] {message}\n")

# Function to run shell commands and log status
def run_command(command, step_name, log_file, critical=False):
    try:
        subprocess.run(command, shell=True, check=True)
        log_status(log_file, f"{step_name}: SUCCESS")
    except subprocess.CalledProcessError:
        log_status(log_file, f"{step_name}: FAILURE")
        print(f"Error: {step_name} failed. Check {log_file} for details.")
        
        if critical:
            exit(1)  # Stop execution if a critical step fails


def extract_contigs(row, fasta_sequences, output_dir):
    ranges = eval(row.iloc[13])
    sequence_id = row.iloc[3]
    sequence = next((seq for seq in fasta_sequences if seq.id == sequence_id), None)
    if sequence is None:
        print(f"Sequence {sequence_id} not found in fasta file.")
        return

    for i, (start, end) in enumerate(ranges):
        exon_name = row['exon_names'][i]
        contig = sequence[start:end]
        contig.id = f"{exon_name}_{sequence_id}_{i+1}"
        contig.description = ""
        file_path = os.path.join(output_dir, f"{exon_name}.fasta")
        with open(file_path, "a") as output_handle:
            SeqIO.write(contig, output_handle, "fasta")

def clean_fasta(row, fasta_sequences, output_dir):
    ranges = eval(row.iloc[13])
    sequence_id = row.iloc[3]
    sequence = next((seq for seq in fasta_sequences if seq.id == sequence_id), None)
    if sequence is None:
        print(f"Sequence {sequence_id} not found in fasta file.")
        return

    for i, (start, end) in enumerate(ranges):
        exon_name = row['exon_names'][i]
        file_path = os.path.join(output_dir, f"{exon_name}.fasta")
        if os.path.exists(file_path):
            os.remove(file_path)

def check_overlap(exon_ranges, start, end, overlap_percentage):
    for (exon_start, exon_end), exon_name in exon_ranges:
        overlap = min(end, exon_end) - max(start, exon_start)
        if overlap > 0 and (overlap / (max(end, exon_end) - min(start, exon_start)) >= overlap_percentage):
            return exon_name
    return None

def process_exon_data(input_dir, gene_name, output_dir, overlap_percentage):
    data_name = os.path.basename(input_dir)
    file_path = os.path.join(input_dir, gene_name, data_name, 'exonerate_stats.tsv')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df = pd.read_csv(file_path, sep='\t')
    end_index = df[df.iloc[:, 0] == 'Hits filtered to remove hits with frameshifts'].index[0]
    df = df.loc[:end_index-1]

    exon_ranges = []
    exon_names = []

    for index, row in df.iterrows():
        ranges = eval(row.iloc[6])
        row_names = []
        for i, (start, end) in enumerate(ranges):
            overlap_name = check_overlap(exon_ranges, start, end, overlap_percentage)
            if overlap_name is None:
                exon_name = f'{data_name}_{gene_name}_exon_{i+1}'
                exon_ranges.append(((start, end), exon_name))
            else:
                exon_name = overlap_name
            row_names.append(exon_name)
        exon_names.append(row_names)
        print(f"Row {index+1}: {len(row_names)} exons and contigs extracted")

    df['exon_names'] = exon_names

    output_file_path = os.path.join(output_dir, f"{data_name}_{gene_name}_exon_split.tsv")
    df.to_csv(output_file_path, sep='\t', index=False)

    fasta_file_path = os.path.join(input_dir, gene_name, f"{gene_name}_contigs.fasta")
    fasta_sequences = list(SeqIO.parse(fasta_file_path, "fasta"))

    for index, row in df.iterrows():
        if row.iloc[9] == "-1":
            row['exon_names'] = row['exon_names'][::-1]
        clean_fasta(row, fasta_sequences, output_dir)

    for index, row in df.iterrows():
        if row.iloc[9] == "-1":
            row['exon_names'] = row['exon_names'][::-1]
        extract_contigs(row, fasta_sequences, output_dir)

# Sequence assembly
def sequence_assembly(threads, read1, read2, mega353, proj_name, log_file, output_hyb):
    # Step 1: Run fastp for trimming
    fastp_cmd = (
        f"fastp -i {read1} -I {read2} -o {read1}.trimmed.fastq.gz -O {read2}.trimmed.fastq.gz "
        f"-j fastp.json -h fastp.html"
    )
    run_command(fastp_cmd, "Sequence Trimming (fastp)", log_file, critical=True)

    # Step 2: Create output directory
    os.makedirs(output_hyb, exist_ok=True)
    log_status(log_file, f"Create Output Directory ({output_hyb}): SUCCESS")

    # Step 3: Run hybpiper for sequence assembly
    hybpiper_cmd = (
        f"hybpiper assemble -t_dna {mega353} "
        f"-r {read1}.trimmed.fastq.gz {read2}.trimmed.fastq.gz "
        f"--prefix {proj_name} --bwa --cpu {threads} -o ./{output_hyb}"
    )
    run_command(hybpiper_cmd, "Sequence Assembly (hybpiper)", log_file, critical=True)


# Exon tree creation
def exon_tree(gene_list, overlapping_rate, proj_name, log_file, output_hyb, output_exon):
    # Step 4: Modify gene list file
    run_command(f"sed s/.fasta//g {gene_list} > gene_list.txt", "Modify Gene List", log_file)

    # Step 5: Create output directory for exon extraction
    os.makedirs(output_exon, exist_ok=True)
    log_status(log_file, f"Create Output Directory ({output_exon}): SUCCESS")
    file_path = os.path.join(output_hyb, proj_name)

    # Read the list of gene names from a file
    with open('gene_list.txt', 'r') as file:
        lines = file.readlines()
    
    lines = [line.strip() for line in lines]
    
    # Step 6: Process the gene data and extract exons
    for gene in lines:  # Assuming lines is a list of gene names
        try:
            # Process exon data for the current gene
            process_exon_data(file_path, gene, output_exon, overlapping_rate)
            # Log the success status
            log_status(log_file, f"Processed Exons for Gene {gene}: SUCCESS")
        except Exception as e:
            # Log the failure status with the exception message
            log_status(log_file, f"Failed to Process Exons for Gene {gene}: {str(e)}: FAILURE")


# Main execution
if __name__ == "__main__":
    # Argument parser for command line inputs
    parser = argparse.ArgumentParser(description="Pipeline for analyzing raw reads from a mixed sample.")
    parser.add_argument("-t", "--threads", required=True, help="Number of threads to use")
    parser.add_argument("-r1", "--read1", required=True, help="Path to the first read file")
    parser.add_argument("-r2", "--read2", required=True, help="Path to the second read file")
    parser.add_argument("-m", "--mega353", default="angiosperms353_v2_interim_targetfile.fasta", help="Path to the mega353 file")
    parser.add_argument("-p", "--proj_name", type=str, required=True, help="Project name")
    parser.add_argument("-g", "--gene_list", required=True, help="Path to the gene list file")
    parser.add_argument("-ov", "--overlapping_rate", type=float, default=0.8, help="Overlapping to consider the same exon (0-1, default 0.8)")
    parser.add_argument("--output_hyb", default="01_hyb_output", help="Custom output folder for hyb results (default is '01_hyb_output')")
    parser.add_argument("--output_exon", default="02_exon_extracted", help="Custom output folder for exon results (default is '02_exon_extracted')")


# Main execution
if __name__ == "__main__":
    # Parse the arguments
    args = parser.parse_args()

    project_name = args.proj_name 

    if not is_valid_project_name(project_name):
        print(f"Error: The project name '{project_name}' is invalid.")
        sys.exit(1)

    # Define the log file name
    log_file = f"{args.proj_name}_01_exon_assembly.log"
    if os.path.exists(log_file):
        os.remove(log_file)

    # Log the input parameters
    log_status(log_file, "Pipeline started with the following parameters:")
    log_status(log_file, f"  Threads: {args.threads}")
    log_status(log_file, f"  Read1: {args.read1}")
    log_status(log_file, f"  Read2: {args.read2}")
    log_status(log_file, f"  Mega353: {args.mega353}")
    log_status(log_file, f"  Project Name: {args.proj_name}")
    log_status(log_file, f"  Gene List: {args.gene_list}")
    log_status(log_file, f"  Output Hyb: {args.output_hyb}")
    log_status(log_file, f"  Output Exon: {args.output_exon}")

    # Run the sequence assembly
    sequence_assembly(args.threads, args.read1, args.read2, args.mega353, args.proj_name, log_file, args.output_hyb)

    # Run the exon tree creation
    exon_tree(args.gene_list, args.overlapping_rate, project_name, log_file, args.output_hyb, args.output_exon)

    log_status(log_file, "Pipeline completed successfully.")
    print(f"Pipeline completed. Check {log_file} for the status of each step.")

