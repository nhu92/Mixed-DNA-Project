import os
import argparse
from Bio import SeqIO

def process_fasta_files(directory):
    valid_files = []
    for filename in os.listdir(directory):
        if filename.endswith(".fasta"):
            file_path = os.path.join(directory, filename)
            all_sequences_valid = True

            for record in SeqIO.parse(file_path, "fasta"):
                if len(record.seq) < 40:
                    all_sequences_valid = False
                    break

            if all_sequences_valid:
                valid_files.append(filename)

    return valid_files

def main():
    parser = argparse.ArgumentParser(description="Process FASTA files")
    parser.add_argument("directory", type=str, help="Directory containing FASTA files")
    parser.add_argument("output_file", type=str, help="Output file to save valid filenames")

    args = parser.parse_args()

    valid_files = process_fasta_files(args.directory)

    with open(args.output_file, 'w') as f:
        for file in valid_files:
            f.write(f"{file}\n")

if __name__ == "__main__":
    main()
