import os
import argparse
import glob

def find_and_list_assembled_genes(source_dir, output_dir):
    # Get a list of all 1st level folder names in the source directory
    folder_names = [folder for folder in os.listdir(source_dir) if os.path.isdir(os.path.join(source_dir, folder))]

    # Initialize a list to store folder names with the specified file
    assembled_gene_folders = []

    for folder_name in folder_names:
        # Check if a file with the pattern "folder name_contigs.fasta" exists
        pattern = os.path.join(source_dir, folder_name, f"{folder_name}_contigs.fasta")
        if glob.glob(pattern):
            assembled_gene_folders.append(folder_name)

    # Write the folder names to the output file
    output_filename = os.path.join(output_dir, f"{os.path.basename(source_dir.rstrip('/'))}_assembled_genes.txt")
    with open(output_filename, 'w') as output_file:
        for folder_name in assembled_gene_folders:
            output_file.write(folder_name + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch assembled genes from HybPiper output")
    parser.add_argument("source_dir", help="Path to the parent of HybPiper output directory")
    parser.add_argument("output_dir", help="Path to the output directory")
    args = parser.parse_args()

    source_dir = args.source_dir
    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    find_and_list_assembled_genes(source_dir, output_dir)
