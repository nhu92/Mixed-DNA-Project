import argparse
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from multiprocessing import Pool

def is_overlapping(args):
    seq1, seq2, min_overlap = args
    aligner = PairwiseAligner()
    alignments = aligner.align(seq1, seq2)
    # Get the number of matching bases in the best alignment
    match_count = sum(a == b for a, b in zip(*alignments[0]))
    return match_count >= min_overlap

def remove_non_overlapping(input_file, output_file, removed_file, min_overlap, num_threads):
    records = list(SeqIO.parse(input_file, "fasta"))
    overlapping_records = []
    removed_records = []
    with Pool(num_threads) as pool:
        for i in range(len(records)):
            args = [(str(records[i].seq), str(records[j].seq), min_overlap) for j in range(i+1, len(records))]
            results = pool.map(is_overlapping, args)
            if any(results):
                overlapping_records.append(records[i])
            else:
                removed_records.append(records[i])
    SeqIO.write(overlapping_records, output_file, "fasta")
    SeqIO.write(removed_records, removed_file, "fasta")
    print(f"Removed {len(removed_records)} sequences")

def main():
    parser = argparse.ArgumentParser(description='Remove non-overlapping sequences from a FASTA file.')
    parser.add_argument('input_file', type=str, help='Input FASTA file')
    parser.add_argument('output_file', type=str, help='Output FASTA file for overlapping sequences')
    parser.add_argument('removed_file', type=str, help='Output FASTA file for removed sequences')
    parser.add_argument('--min_overlap', type=int, default=100,
                        help='Minimum number of matching bases for two sequences to be considered overlapping')
    parser.add_argument('--num_threads', type=int, default=1,
                        help='Number of threads to use')
    args = parser.parse_args()
    remove_non_overlapping(args.input_file, args.output_file, args.removed_file,
                           args.min_overlap, args.num_threads)

if __name__ == "__main__":
    main()
