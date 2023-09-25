import os
import argparse
from Bio import SeqIO

def filter_gene_by_number(input_dir, gene_list, max_sequences, output_dir):
    for gene_name in gene_list:
        fasta_file = os.path.join(input_dir, f"{gene_name}_reduced.fasta")
        output_file = os.path.join(output_dir, f"{gene_name}_filtered.fasta")
        
        if not os.path.exists(fasta_file):
            print(f"Warning: {fasta_file} does not exist.")
            continue
        
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        num_sequences = len(sequences)
        
        if num_sequences <= max_sequences:
            SeqIO.write(sequences, output_file, "fasta")
        else:
            print(f"Removed {gene_name}_reduced.fasta with {num_sequences} sequences.")
            continue

def main():
    parser = argparse.ArgumentParser(description="Filter FASTA files based on the number of sequences.")
    parser.add_argument("gene_list", help="Path to the gene name list file")
    parser.add_argument("input_dir", help="Path to the input directory containing FASTA files")
    parser.add_argument("max_sequences", type=int, help="Maximum number of sequences to keep in a file")
    parser.add_argument("output_dir", help="Path to the output directory for filtered files")
    
    args = parser.parse_args()
    
    with open(args.gene_list, "r") as gene_list_file:
        gene_list = [line.strip() for line in gene_list_file]
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    filter_gene_by_number(args.input_dir, gene_list, args.max_sequences, args.output_dir)

if __name__ == "__main__":
    main()
