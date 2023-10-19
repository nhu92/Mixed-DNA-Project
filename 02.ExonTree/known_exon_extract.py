import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

def extract_exons(input_fasta, input_gff, gene_name, species_name, output_directory):
    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    gene_sequences = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
    gff_data = open(input_gff, "r")

    exon_number = 1

    for line in gff_data:
        if not line.strip() or line.startswith("#"):
            continue
        columns = line.strip().split("\t")
        feature_type = columns[2]

        if feature_type == "exon":
            start = int(columns[3]) - 1  # Convert to 0-based index
            end = int(columns[4])
            seqname = columns[0]
            gene_sequence = gene_sequences.get(seqname)

            if gene_sequence is not None:
                exon_sequence = gene_sequence.seq[start:end]

                # Create a SeqRecord for the exon
                exon_record = SeqRecord(exon_sequence, id=f"{species_name}.{gene_name}.exon{exon_number}", description="")
                exon_number += 1

                # Save the exon as a FASTA file
                exon_filename = os.path.join(output_directory, f"{species_name}.{gene_name}.exon{exon_number - 1}.fasta")
                SeqIO.write(exon_record, exon_filename, "fasta")

    gff_data.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract exons from a GFF file and save them as FASTA files.")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("input_gff", help="Input GFF file")
    parser.add_argument("gene_name", help="Gene name")
    parser.add_argument("species_name", help="Species name")
    parser.add_argument("output_directory", help="Output directory")

    args = parser.parse_args()

    extract_exons(args.input_fasta, args.input_gff, args.gene_name, args.species_name, args.output_directory)
