import os
import pandas as pd
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--input_dir', required=True)
    parser.add_argument('--gene_name', required=True)
    parser.add_argument('--output_dir', required=True)
    return parser.parse_args()

def create_output_dir(output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

def get_exonerate_stats_file(args):
    data_name = os.path.basename(args.input_dir)
    return os.path.join(args.input_dir, args.gene_name, data_name, 'exonerate_stats.tsv')

def get_fasta_file(args):
    return os.path.join(args.input_dir, args.gene_name, f"{args.gene_name}_contigs.fasta")

def get_df(exonerate_stats_file):
    df = pd.read_csv(exonerate_stats_file, sep='\t')
    return df.loc[:df[df.iloc[:,0] == "Hits with subsumed hits removed"].index[0] - 1]

def write_sequences(df, fasta_file, args):
    data_name = os.path.basename(args.input_dir)
    output_file = os.path.join(args.output_dir, f"{data_name}_{args.gene_name}.fasta")

    with open(output_file, 'w') as f:
        for index, row in df.iterrows():
            exon_ranges_str = row.iloc[6]
            sequence_name = row.iloc[3]
            exon_ranges = eval(exon_ranges_str)

            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id == sequence_name:
                    sequence_parts = []
                    for start, end in exon_ranges:
                        sequence_parts.append(str(record.seq[start:end]))
                    sequence = ("-"*10).join(sequence_parts)
                    f.write(f">{data_name}_{args.gene_name}_{sequence_name}\n{sequence}\n")

def main():
    args = get_args()

    create_output_dir(args.output_dir)

    exonerate_stats_file = get_exonerate_stats_file(args)

    fasta_file = get_fasta_file(args)

    df = get_df(exonerate_stats_file)

    write_sequences(df, fasta_file, args)

if __name__ == "__main__":
    main()
