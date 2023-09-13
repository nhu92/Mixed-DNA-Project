import argparse
from Bio import SeqIO
from collections import Counter
from scipy.stats import chisquare
import numpy as np

def compute_frequencies(sequences):
    counter = Counter()
    total = 0
    for sequence in sequences:
        sequence = [base for base in sequence if base in 'ATCG']
        counter.update(sequence)
        total += len(sequence)
    return {base: count / total for base, count in counter.items()}

def composition_test(sequence, expected_frequencies):
    observed_frequencies = compute_frequencies([sequence])
    observed = np.array([observed_frequencies.get(base, 0) * len(sequence) for base in expected_frequencies])
    expected = np.array([expected_frequencies[base] * len(sequence) for base in expected_frequencies])
    chi2, p = chisquare(observed, expected)
    return p > 0.05, p

def main(input_file, output_file):
    with open(input_file, 'r') as file:
        sequences = [record for record in SeqIO.parse(file, 'fasta')]
        expected_frequencies = compute_frequencies([str(record.seq) for record in sequences])
        passed_sequences = []
        for record in sequences:
            passed, p_value = composition_test(str(record.seq), expected_frequencies)
            result = 'Passed' if passed else 'Failed'
            print(f'Sequence {record.id}: {result}, p-value: {p_value}')
            if passed:
                passed_sequences.append(record)
        SeqIO.write(passed_sequences, output_file, 'fasta')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform a composition test on fasta sequences')
    parser.add_argument('input_file', type=str, help='Input fasta file')
    parser.add_argument('output_file', type=str, help='Output fasta file')

    args = parser.parse_args()
    main(args.input_file, args.output_file)
