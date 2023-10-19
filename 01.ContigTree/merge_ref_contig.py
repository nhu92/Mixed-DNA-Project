import os
import argparse

def merge_gene_fasta(fna_dir, fasta_dir, output_dir, gene_list_file):
    with open(gene_list_file, "r") as gene_list:
        gene_names = [line.strip() for line in gene_list]

    for gene_name in gene_names:
        fna_files = []
        fasta_files = []

        # Search for FNA files in the FNA directory (case-insensitive match)
        for filename in os.listdir(fna_dir):
            if filename.lower().startswith(gene_name.lower()) and filename.lower().endswith(".fna"):
                fna_files.append(os.path.join(fna_dir, filename))

        # Search for FASTA files in the FASTA directory (case-insensitive match)
        for filename in os.listdir(fasta_dir):
            if filename.lower().startswith(gene_name.lower()) and filename.lower().endswith(".fasta"):
                fasta_files.append(os.path.join(fasta_dir, filename))

        if fna_files and fasta_files:
            merged_sequences = []

            # Read and merge all matching FNA files
            for fna_file in fna_files:
                with open(fna_file, "r") as fna:
                    merged_sequences.extend(fna.readlines())

            # Read and merge all matching FASTA files
            for fasta_file in fasta_files:
                with open(fasta_file, "r") as fasta:
                    merged_sequences.extend(fasta.readlines())

            # Write merged sequences to the output file
            output_file = os.path.join(output_dir, f"{gene_name}_ref_contig_merged.fasta")
            with open(output_file, "w") as merged_file:
                merged_file.writelines(merged_sequences)

            print(f"Merged {len(fna_files) + len(fasta_files)} files for {gene_name}.")

        else:
            print(f"Files for gene {gene_name} not found in the input directories.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge gene sequences from FNA and FASTA files.")
    parser.add_argument("fna_dir", help="Input directory containing FNA files.")
    parser.add_argument("fasta_dir", help="Input directory containing FASTA files.")
    parser.add_argument("output_dir", help="Output directory for merged files.")
    parser.add_argument("gene_list_file", help="File containing a list of gene names, one per line.")

    args = parser.parse_args()
    fna_dir = args.fna_dir
    fasta_dir = args.fasta_dir
    output_dir = args.output_dir
    gene_list_file = args.gene_list_file

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    merge_gene_fasta(fna_dir, fasta_dir, output_dir, gene_list_file)
