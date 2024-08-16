import subprocess
import os
import argparse

# Function to run shell commands and log output
def run_command(command, log_file):
    try:
        with open(log_file, 'a') as log:
            log.write(f"Running command: {command}\n")
            subprocess.run(command, shell=True, check=True, stdout=log, stderr=log)
            log.write("Command finished successfully.\n\n")
    except subprocess.CalledProcessError as e:
        with open(log_file, 'a') as log:
            log.write(f"Error running command: {command}\n")
            log.write(f"{str(e)}\n\n")
        print(f"Error: Command failed. Check {log_file} for details.")

# Sequence assembly
def sequence_assembly(threads, read1, read2, mega353, proj_name, log_file):
    # Step 1: Run fastp for trimming
    fastp_cmd = (
        f"fastp -i {read1} -I {read2} -o {read1}.trimmed.fastq.gz -O {read2}.trimmed.fastq.gz "
        f"-j fastp.json -h fastp.html"
    )
    run_command(fastp_cmd, log_file)

    # Step 2: Create output directory
    os.makedirs("hyb_output", exist_ok=True)

    # Step 3: Run hybpiper for sequence assembly
    hybpiper_cmd = (
        f"hybpiper assemble -t_dna {mega353} "
        f"-r {read1}.trimmed.fastq.gz {read2}.trimmed.fastq.gz "
        f"--prefix {proj_name} --bwa --cpu {threads} -o ./hyb_output"
    )
    run_command(hybiper_cmd, log_file)

# Exon tree creation
def exon_tree(gene_list, proj_name, log_file):
    # Step 4: Modify gene list file
    run_command(f"sed s/.fasta//g {gene_list} > gene_list.txt", log_file)

    # Step 5: Create output directory for exon extraction
    os.makedirs("./exon_extracted", exist_ok=True)

    # Step 6: Extract exons for each gene in the list
    with open("gene_list.txt", "r") as file:
        for line in file:
            line = line.strip()
            split_exon_cmd = (
                f"python split_exon_extract.py ./hyb_output/{proj_name} {line} ./exon_extracted 0.8"
            )
            run_command(split_exon_cmd, log_file)

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

    # Run the sequence assembly
    sequence_assembly(args.threads, args.read1, args.read2, args.mega353, args.proj_name, log_file)

    # Run the exon tree creation
    exon_tree(args.gene_list, args.proj_name, log_file)
