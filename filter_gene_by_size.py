import os
import argparse

def get_gene_file_sizes(input_directory):
    gene_file_sizes = {}
    for filename in os.listdir(input_directory):
        if filename.endswith("_ref_contig_merged.fasta"):
            full_path = os.path.join(input_directory, filename)
            file_size_kb = os.path.getsize(full_path) / 1024  # Convert to KB
            gene_name = filename.replace("_ref_contig_merged.fasta", "")
            gene_file_sizes[gene_name] = file_size_kb
    return gene_file_sizes

def filter_genes_by_size(gene_list, max_size_kb, input_directory):
    gene_file_sizes = get_gene_file_sizes(input_directory)
    removed_genes = []
    filtered_genes = []
    for gene in gene_list:
        if gene in gene_file_sizes and gene_file_sizes[gene] > max_size_kb:
            removed_genes.append((gene, gene_file_sizes[gene]))
        else:
            filtered_genes.append(gene)
    return removed_genes, filtered_genes

def main():
    parser = argparse.ArgumentParser(description="Filter gene names by file size.")
    parser.add_argument("gene_list_file", help="Path to the gene name list file")
    parser.add_argument("input_directory", help="Directory containing gene files")
    parser.add_argument("max_size_kb", type=float, help="Maximum file size in KB")
    parser.add_argument("output_file", help="Output file to store the filtered gene names")
    
    args = parser.parse_args()

    with open(args.gene_list_file, 'r') as gene_list_file:
        gene_list = [line.strip() for line in gene_list_file.readlines()]

    removed_genes, filtered_genes = filter_genes_by_size(gene_list, args.max_size_kb, args.input_directory)

    with open(args.output_file, 'w') as output_file:
        for gene in filtered_genes:
            output_file.write(f"{gene}\n")

    for gene, size in removed_genes:
        print(f"Removed Gene: {gene}, File Size: {size} KB")

if __name__ == "__main__":
    main()
