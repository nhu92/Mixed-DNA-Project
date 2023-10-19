import os
import argparse
from Bio import SeqIO

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter FASTA files based on the number of sequences")
    parser.add_argument("gene_list", help="Path to the gene name list file")
    parser.add_argument("input_directory", help="Path to the input directory containing FASTA files")
    parser.add_argument("output_directory", help="Path to the output directory for unremoved gene list")
    parser.add_argument("max_sequences", type=int, help="Maximum number of sequences allowed in a FASTA file")
    return parser.parse_args()

def filter_fasta_files(gene_list_path, input_dir, output_dir, max_sequences):
    unremoved_genes = []

    # Read the gene name list
    with open(gene_list_path, 'r') as gene_list_file:
        gene_names = set(line.strip() for line in gene_list_file)

    # Process each FASTA file in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith("_reduced.fasta"):
            fasta_path = os.path.join(input_dir, filename)
            gene_name = filename.replace("_reduced.fasta", "")

            if gene_name not in gene_names:
                continue

            # Count sequences in the FASTA file
            num_sequences = sum(1 for record in SeqIO.parse(fasta_path, "fasta"))

            if num_sequences > max_sequences:
                print(f"Removing {gene_name} as it contains {num_sequences} sequences (more than {max_sequences})")
            else:
                unremoved_genes.append(gene_name)

    # Write unremoved gene list to the output directory
    with open(output_dir, 'w') as output_file:
        output_file.write("\n".join(unremoved_genes))

    print(f"Unremoved genes written to {output_dir}")

if __name__ == "__main__":
    args = parse_arguments()
    filter_fasta_files(args.gene_list, args.input_directory, args.output_directory, args.max_sequences)
