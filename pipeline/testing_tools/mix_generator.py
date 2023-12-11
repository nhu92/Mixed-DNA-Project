import argparse
import numpy as np
import os

def generate_poisson_number(mean=5, min_value=2, max_value=12):
    while True:
        number = np.random.poisson(mean)
        if min_value <= number <= max_value:
            return number

def sample_species_numbers(n_mixed, total_species=30):
    return np.random.choice(range(1, total_species + 1), n_mixed, replace=False)

def merge_files(species_numbers, input_dir, output_dir, read_type):
    files_to_merge = [f"{str(species).zfill(2)}.{read_type}.fastq.gz" for species in species_numbers]
    merged_filename = "x".join([str(species).zfill(2) for species in species_numbers]) + f".{read_type}.fastq.gz"
    with open(os.path.join(output_dir, merged_filename), 'wb') as merged_file:
        for file in files_to_merge:
            with open(os.path.join(input_dir, file), 'rb') as f:
                merged_file.write(f.read())

def main():
    parser = argparse.ArgumentParser(description="Merge species read files.")
    parser.add_argument("input_dir", type=str, help="Directory containing read files.")
    parser.add_argument("output_dir", type=str, help="Directory to save merged read files.")
    args = parser.parse_args()

    n_mixed = generate_poisson_number()
    species_numbers = sample_species_numbers(n_mixed)

    merge_files(species_numbers, args.input_dir, args.output_dir, "R1")
    merge_files(species_numbers, args.input_dir, args.output_dir, "R2")

    print(f"Merged files for species: {'x'.join([str(species).zfill(2) for species in species_numbers])}")

if __name__ == "__main__":
    main()

