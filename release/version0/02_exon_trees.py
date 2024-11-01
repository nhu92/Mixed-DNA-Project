import subprocess
import os
import argparse
from datetime import datetime

# Function to log status with timestamp and flush immediately
def log_status(log_file, message):
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(log_file, 'a') as log:
        log.write(f"[{timestamp}] {message}\n")
        log.flush()

def run_command(command, step_name, log_file, critical=False):
    try:
        subprocess.run(command, shell=True, check=True)
        log_status(log_file, f"{step_name}: SUCCESS")
    except subprocess.CalledProcessError:
        log_status(log_file, f"{step_name}: FAILURE")
        print(f"Error: {step_name} failed. Check {log_file} for details.")
        
        if critical:
            exit(1)

def get_first_sequence_length(file_path):
    """Reads the first sequence from a FASTA file and returns its length."""
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith('>'):  # Skip header lines
                return len(line.strip())  # Return length of first sequence line

# Function to perform the batch run to generate trees for all genes
def generate_trees(threads, input_exon, ref_alignment, gene_list, proj_name, output_dir, log_file, exon_min_size=80):
    os.makedirs(output_dir, exist_ok=True)
    log_status(log_file, f"Created directory {output_dir}")

    with open(gene_list, 'r') as genes:
        for gene_name_shorter in genes:
            gene_name_shorter = gene_name_shorter.strip()
            
            # Non-critical ls command: failure won't stop the pipeline
            run_command(f"ls {input_exon}/*{gene_name_shorter}*.fasta > ./exon.list.txt", 
                        f"List exons for {gene_name_shorter}", log_file, critical=False)

            # Check if the exon list file is empty
            if not os.path.exists('./exon.list.txt') or os.path.getsize('./exon.list.txt') == 0:
                log_status(log_file, f"No exon files found for {gene_name_shorter}, skipping.")
                if os.path.exists('./exon.list.txt'):
                    os.remove('./exon.list.txt')
                continue

            with open('./exon.list.txt', 'r') as exons:
                i = 1
                for exon_file in exons:
                    exon_file = exon_file.strip()

                    # Get the length of the first sequence in the exon file
                    exon_length = get_first_sequence_length(exon_file)
                    if exon_length < exon_min_size:
                        log_status(log_file, f"Skipping {exon_file} (sequence length: {exon_length} < {exon_min_size})")
                        continue

                    # MAFFT alignment and trimming
                    mafft_cmd = (
                        f"mafft --preservecase --maxiterate 1000 --localpair --adjustdirection "
                        f"--thread {threads} --addfragments {exon_file} {ref_alignment}/{gene_name_shorter}.fasta "
                        f"> {output_dir}/{gene_name_shorter}_exon_{i}_aligned.fasta"
                    )
                    run_command(mafft_cmd, f"MAFFT alignment for {gene_name_shorter} exon {i}", log_file)

                    trimal_cmd = (
                        f"trimal -in {output_dir}/{gene_name_shorter}_exon_{i}_aligned.fasta "
                        f"-out {output_dir}/{gene_name_shorter}_exon_{i}_trimmed.fasta -gt 0.5"
                    )
                    run_command(trimal_cmd, f"Trim alignment for {gene_name_shorter} exon {i}", log_file)

                    # Tree construction
                    fasttree_cmd = (
                        f"fasttree -gtr -gamma -nt {output_dir}/{gene_name_shorter}_exon_{i}_trimmed.fasta "
                        f"> {output_dir}/{gene_name_shorter}_exon_{i}.tre"
                    )
                    run_command(fasttree_cmd, f"Tree construction for {gene_name_shorter} exon {i}", log_file)

                    i += 1

            os.remove('./exon.list.txt')
            log_status(log_file, "Removed temporary file exon.list.txt")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline for generating phylogenetic trees from exon data.")
    parser.add_argument("-t", "--threads", required=True, help="Number of threads to use")
    parser.add_argument("-e", "--input_exon", type=str, default="02_exon_extracted", help="Directory of extracted exon sequences")
    parser.add_argument("-r", "--ref_alignment", type=str, default="ref", help="Directory of reference alignments")
    parser.add_argument("-g", "--gene_list", type=str, default="gene_list.txt", help="Path to the gene list file")
    parser.add_argument("-p", "--proj_name", required=True, help="Project name")
    parser.add_argument("-o", "--output_dir", type=str, default="03_phylo_results", help="Directory for output phylogenetic trees")
    parser.add_argument("-m", "--exon_min_size", type=int, default=80, help="Minimum exon size to include in analysis")

    args = parser.parse_args()

    log_file = f"{args.proj_name}_02_exons_phylo.log"
    if os.path.exists(log_file):
        os.remove(log_file)

    log_status(log_file, "Pipeline started with the following parameters:")
    log_status(log_file, f"  Threads: {args.threads}")
    log_status(log_file, f"  Input Exon Directory: {args.input_exon}")
    log_status(log_file, f"  Reference Alignment Directory: {args.ref_alignment}")
    log_status(log_file, f"  Gene List: {args.gene_list}")
    log_status(log_file, f"  Project Name: {args.proj_name}")
    log_status(log_file, f"  Output Directory: {args.output_dir}")
    log_status(log_file, f"  Exon Minimum Size: {args.exon_min_size}")

    generate_trees(args.threads, args.input_exon, args.ref_alignment, args.gene_list, args.proj_name, args.output_dir, log_file, args.exon_min_size)

    log_status(log_file, "Pipeline completed successfully.")
    print(f"Pipeline completed. Check {log_file} for the status of each step.")
