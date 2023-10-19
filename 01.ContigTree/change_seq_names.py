import os
import shutil
import argparse
from Bio import SeqIO

# Function to generate a list of folder names in the 1st level of the source directory
def generate_folder_list(source_directory):
    folder_list = [d for d in os.listdir(source_directory) if os.path.isdir(os.path.join(source_directory, d))]
    return folder_list

# Function to copy and modify FASTA files
def copy_and_modify_fasta(folder_name, input_directory, source_directory):
    source_dir_name = os.path.basename(source_directory.rstrip('/'))
    source_fasta = os.path.join(source_directory, folder_name, f"{folder_name}_contigs.fasta")

    # Check if the source file exists
    if not os.path.exists(source_fasta):
        print(f"Skipping {folder_name}: Source file does not exist.")
        return

    output_filename = f"{folder_name}_{source_dir_name}_contigs.fasta"
    output_fasta = os.path.join(input_directory, output_filename)

    order_number = 1
    records = []

    with open(source_fasta, "r") as source_file:
        for record in SeqIO.parse(source_file, "fasta"):
            new_id = f"{folder_name}_{source_dir_name}_{order_number}"
            record.id = new_id
            record.description = ""
            records.append(record)
            order_number += 1

    with open(output_fasta, "w") as output_file:
        SeqIO.write(records, output_file, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Copy and modify names of contigs files from HybPiper.")
    parser.add_argument("source_directory", help="Path to the parent directory of HybPiper results containing subdirectories and FASTA files.")
    parser.add_argument("input_directory", help="Path to the output directory where modified FASTA files will be copied.")
    args = parser.parse_args()

    source_directory = args.source_directory
    input_directory = args.input_directory
    source_species = os.path.basename(source_directory.rstrip('/'))
    print(f"Currently modifying: {source_species}.")
    # Generate a list of folder names in the 1st level of the source directory
    folders = generate_folder_list(source_directory)

    # Generate an output directroy
    if not os.path.exists(input_directory):
        os.makedirs(input_directory)
        
    # Iterate over folders, copy and modify FASTA files
    for folder in folders:
        copy_and_modify_fasta(folder, input_directory, source_directory)

    print("FASTA files copied and modified successfully.")
