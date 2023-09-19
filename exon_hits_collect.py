import os
import shutil
import argparse
import pandas as pd
from Bio import SeqIO

# Function to generate a list of folder names in the 1st level of the source directory
def generate_folder_list(source_directory):
    folder_list = [d for d in os.listdir(source_directory) if os.path.isdir(os.path.join(source_directory, d))]
    return folder_list

# Function to process TSV file and return sequence info
def process_tsv(tsv_file):
    # Check if the source file exists
    if not os.path.exists(tsv_file):
        print(f"Skipping {tsv_file}: Source file does not exist.")
        exit()

    df = pd.read_csv(tsv_file, sep='\t')

    index = df[df.iloc[:, 0] == 'Hits with subsumed hits removed'].index[0]
    df = df.iloc[:index]

    sequence_info = []
    for index, row in df.iterrows():
        sequence_name = row.iloc[3]  # Fourth column
        start, end = map(int, row.iloc[11].strip('()').split(', '))  # Twelveth column
        sequence_info.append((sequence_name, start, end))

    return sequence_info

def trim_fasta(input_fasta, output_fasta, sequence_info, folder_name, source_dir_name):
    order_number = 1

    with open(output_fasta, 'w') as output_file:
        for record in SeqIO.parse(input_fasta, "fasta"):
            sequence_name = record.id
            sequence = str(record.seq)

            for name, start, end in sequence_info:
                if name == sequence_name:
                    trimmed_sequence = sequence[start:end]
                    new_id = f"{folder_name}_{source_dir_name}_{order_number}"
                    output_file.write(f'>{new_id}\n') # testing
                    output_file.write(f'{trimmed_sequence}\n')
                    order_number += 1

# Function to copy and modify FASTA files
def copy_and_modify_fasta(folder_name, input_directory, source_directory, sequence_info):
    source_dir_name = os.path.basename(source_directory.rstrip('/'))
    source_fasta = os.path.join(source_directory, folder_name, f"{folder_name}_contigs.fasta")

    # Check if the source file exists
    if not os.path.exists(source_fasta):
        print(f"Skipping {folder_name}: Source file does not exist.")
        return

    output_filename = f"{folder_name}_{source_dir_name}_contigs.fasta"
    output_fasta = os.path.join(input_directory, output_filename)

    trim_fasta(source_fasta, output_fasta, sequence_info, folder_name, source_dir_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Copy, trim, and modify names of contigs files from HybPiper.")
    parser.add_argument("source_directory", help="Path to the parent directory of HybPiper results containing subdirectories and FASTA files.")
    parser.add_argument("input_directory", help="Path to the output directory where modified FASTA files will be copied.")
    args = parser.parse_args()

    source_directory = args.source_directory
    input_directory = args.input_directory

    # Generate a list of folder names in the 1st level of the source directory
    folders = generate_folder_list(source_directory)

    # Generate an output directroy
    if not os.path.exists(input_directory):
        os.makedirs(input_directory)

    # Iterate over folders, copy, trim, and modify FASTA files
    for folder in folders:
        tsv_file = os.path.join(source_directory, folder, os.path.basename(source_directory.rstrip('/')), "exonerate_stats.tsv")
        sequence_info = process_tsv(tsv_file)
        copy_and_modify_fasta(folder, input_directory, source_directory, sequence_info)

    print("FASTA files copied, modified, and trimmed successfully.")
