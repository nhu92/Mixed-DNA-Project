import argparse
from Bio import SeqIO

def extract_sequence_names(input_fasta):
    sequence_names = set()
    with open(input_fasta, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Remove '_R_' tags from the beginning of the sequence name
            sequence_name = record.id.lstrip('_R_')
            sequence_names.add(sequence_name)
    return sequence_names

def extract_sequences(input_fasta, sequence_names, output_fasta):
    with open(input_fasta, "r") as input_file, open(output_fasta, "w") as output_file:
        for record in SeqIO.parse(input_file, "fasta"):
            sequence_name = record.id.lstrip('_R_')
            if sequence_name in sequence_names:
                SeqIO.write(record, output_file, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Extract and filter sequences from FASTA files.")
    parser.add_argument("input_fasta_names", help="Input FASTA file containing sequence names")
    parser.add_argument("input_fasta_sequences", help="Input FASTA file containing sequences")
    parser.add_argument("output_fasta", help="Output FASTA file")
    
    args = parser.parse_args()

    sequence_names = extract_sequence_names(args.input_fasta_names)
    extract_sequences(args.input_fasta_sequences, sequence_names, args.output_fasta)

if __name__ == "__main__":
    main()
