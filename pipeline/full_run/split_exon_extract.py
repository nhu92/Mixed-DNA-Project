import argparse
import os
import pandas as pd
from intervaltree import Interval, IntervalTree
from Bio import SeqIO

def extract_contigs(row, fasta_sequences, output_dir):
    ranges = eval(row.iloc[13])
    sequence_id = row.iloc[3]
    sequence = next((seq for seq in fasta_sequences if seq.id == sequence_id), None)
    if sequence is None:
        print(f"Sequence {sequence_id} not found in fasta file.")
        return

    for i, (start, end) in enumerate(ranges):
        exon_name = row['exon_names'][i]
        contig = sequence[start:end]
        # Include both the exon name and a unique identifier in the contig ID
        contig.id = f"{exon_name}_{sequence_id}_{i+1}"
        contig.description = ""

        # Write to a fasta file
        file_path = os.path.join(output_dir, f"{exon_name}.fasta")
        with open(file_path, "a") as output_handle:
            SeqIO.write(contig, output_handle, "fasta")


def clean_fasta(row, fasta_sequences, output_dir):
    ranges = eval(row.iloc[13])
    sequence_id = row.iloc[3]
    sequence = next((seq for seq in fasta_sequences if seq.id == sequence_id), None)
    if sequence is None:
        print(f"Sequence {sequence_id} not found in fasta file.")
        return

    for i, (start, end) in enumerate(ranges):
        exon_name = row['exon_names'][i]
        file_path = os.path.join(output_dir, f"{exon_name}.fasta")
        if os.path.exists(file_path):
            os.remove(file_path)


def check_overlap(exon_ranges, start, end, overlap_percentage):
    for (exon_start, exon_end), exon_name in exon_ranges:
        overlap = min(end, exon_end) - max(start, exon_start)
        if overlap > 0 and (overlap / (max(end, exon_end) - min(start, exon_start)) >= overlap_percentage):
            return exon_name
    return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_dir')
    parser.add_argument('gene_name')
    parser.add_argument('output_dir')
    parser.add_argument('overlap_percentage', type=float)
    args = parser.parse_args()

    data_name = os.path.basename(args.input_dir)
    file_path = os.path.join(args.input_dir, args.gene_name, data_name, 'exonerate_stats.tsv')

    # Make the output directory
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    df = pd.read_csv(file_path, sep='\t')
    # Find the index of the row where the first cell is "Hits with subsumed hits removed"
    end_index = df[df.iloc[:, 0] == 'Hits with subsumed hits removed'].index[0]

    # Only keep the rows above this index
    df = df.loc[:end_index-1]

    exon_ranges = []
    exon_names = []

    for index, row in df.iterrows():
        ranges = eval(row.iloc[6])
        row_names = []
        for i, (start, end) in enumerate(ranges):
            overlap_name = check_overlap(exon_ranges, start, end, args.overlap_percentage)
            if overlap_name is None:
                exon_name = f'{data_name}_{args.gene_name}_exon_{i+1}'
                exon_ranges.append(((start, end), exon_name))
            else:
                exon_name = overlap_name
            row_names.append(exon_name)
        exon_names.append(row_names)

        # Print the number of exons and contigs extracted for each exon
        print(f"Row {index+1}: {len(row_names)} exons and contigs extracted")

    df['exon_names'] = exon_names

    output_file_path = os.path.join(args.output_dir, f"{data_name}_{args.gene_name}_exon_split.tsv")
    df.to_csv(output_file_path, sep='\t', index=False)

    fasta_file_path = os.path.join(args.input_dir, args.gene_name, f"{args.gene_name}_contigs.fasta")
    fasta_sequences = list(SeqIO.parse(fasta_file_path, "fasta"))

    for index, row in df.iterrows():
        if row.iloc[9] == "-1":
            row['exon_names'] = row['exon_names'][::-1]
        clean_fasta(row, fasta_sequences, args.output_dir)


    for index, row in df.iterrows():
        if row.iloc[9] == "-1":
            row['exon_names'] = row['exon_names'][::-1]
        extract_contigs(row, fasta_sequences, args.output_dir)
if __name__ == '__main__':
    main()

