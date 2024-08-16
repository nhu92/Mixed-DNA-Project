import subprocess
import os
import argparse
from datetime import datetime

# Function to log status with timestamp
def log_status(log_file, message):
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(log_file, 'a') as log:
        log.write(f"[{timestamp}] {message}\n")

# Function to run shell commands and log status
def run_command(command, step_name, log_file):
    try:
        subprocess.run(command, shell=True, check=True)
        log_status(log_file, f"{step_name}: SUCCESS")
    except subprocess.CalledProcessError:
        log_status(log_file, f"{step_name}: FAILURE")
        print(f"Error: {step_name} failed. Check {log_file} for details.")
        exit(1)  # Stop execution if a step fails

# Sequence assembly
def sequence_assembly(threads, read1, read2, mega353, proj_name, log_file):
    # Step 1: Run fastp for trimming
    fastp_cmd = (
        f"fastp -i {read1} -I {read2} -o {read1}.trimmed.fastq.gz -O {read2}.trimmed.fastq.gz "
        f"-j fastp.json -h fastp.html"
    )
    run_command(fastp_cmd, "Sequence Trimming (fastp)", log_file)

    # Step 2: Create output directory
    os.makedirs("01_hyb_output", exist_ok=True)
    log_status(log_file, "Create Output Directory (hyb_output): SUCCESS")

    # Step 3: Run hybpiper for sequence assembly
    hybpiper_cmd = (
        f"hybpiper assemble -t_dna {mega353} "
        f"-r {read1}.trimmed.fastq.gz {read2}.trimmed.fastq.gz "
        f"--prefix {proj_name} --bwa --cpu {threads} -o ./01_hyb_output"
    )
    run_command(hybpiper_cmd, "Sequence Assembly (hybpiper)", log_file)

# Exon tree creation
def exon_tree(gene_list, proj_name, log_file):
    # Step 4: Modify gene list file
    run_command(f"sed s/.fasta//g {gene_list} > gene_list.txt", "Modify Gene List", log_file)

    # Step 5: Create output directory for exon extraction
    os.makedirs("./02_exon_extracted", exist_ok=True)
    log_status(log_file, "Create Output Directory (exon_extracted): SUCCESS")

    # Step 6: Extract exons for each gene in the list
    with open("gene_list.txt", "r") as file:
        for line in file:
            line = line.strip()
            split_exon_cmd = (
                f"python split_exon_extract.py ./01_hyb_output/{proj_name} {line} ./02_exon_extracted 0.8"
            )
            run_command(split_exon_cmd, f"Extract Exons for Gene {line}", log_file)

# Main execution
if __name__ == "__main__":
    # Argument parser for command line inputs
    parser = argparse.ArgumentParser(description="Pipeline for analyzing raw reads from a mixed sample.")
    parser.add_argument("-t", "--threads", required=True, help="Number of threads to use")
    parser.add_argument("-r1", "--read1", required=True, help="Path to the first read file")
    parser.add_argument("-r2", "--read2", required=True, help="Path to the second read file")
    parser.add_argument("-m", "--mega353", default="angiosperms353_v2_interim_targetfile.fasta", help="Path to the mega353 file")
    parser.add_argument("-p", "--proj_name", required=True, help="Project name")
    parser.add_argument("-g", "--gene_list", required=True, help="Path to the gene list file")

    # Parse the arguments
    args = parser.parse_args()

    # Define the log file name
    log_file = f"{args.proj_name}_01_exon_assembly.out"

    # Log the input parameters
    log_status(log_file, "Pipeline started with the following parameters:")
    log_status(log_file, f"  Threads: {args.threads}")
    log_status(log_file, f"  Read1: {args.read1}")
    log_status(log_file, f"  Read2: {args.read2}")
    log_status(log_file, f"  Mega353: {args.mega353}")
    log_status(log_file, f"  Project Name: {args.proj_name}")
    log_status(log_file, f"  Gene List: {args.gene_list}")

    # Run the sequence assembly
    sequence_assembly(args.threads, args.read1, args.read2, args.mega353, args.proj_name, log_file)

    # Run the exon tree creation
    exon_tree(args.gene_list, args.proj_name, log_file)

    log_status(log_file, "Pipeline completed successfully.")
    print(f"Pipeline completed. Check {log_file} for the status of each step.")
