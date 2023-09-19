import os
import argparse

# Function to find shared genes among species
def find_shared_genes(input_dir, species_file, output_dir):
    # Read species names from the species file
    with open(species_file, "r") as species_file:
        species_names = [line.strip() for line in species_file]

    # Initialize a dictionary to store genes for each species
    species_genes = {species: set() for species in species_names}

    # Initialize a dictionary to store total gene count for each species
    total_gene_count = {species: 0 for species in species_names}

    # List files in the input directory
    files = os.listdir(input_dir)

    # Iterate through files and extract genes
    for file_name in files:
        if file_name.endswith("_assembled_genes.txt"):
            species_name = file_name.split("_")[0]
            file_path = os.path.join(input_dir, file_name)

            with open(file_path, "r") as file:
                gene_set = species_genes[species_name]
                for line in file:
                    gene_name = line.strip()
                    gene_set.add(gene_name)
                    total_gene_count[species_name] += 1

    # Find the shared genes among all species
    shared_genes = set.intersection(*species_genes.values())

    # Print gene statistics for each species
    for species_name in species_names:
        not_shared_count = total_gene_count[species_name] - len(shared_genes)
        print(f"{species_name}: {total_gene_count[species_name]} genes assembled, {not_shared_count} is/are unshared and removed")

    # Write the shared genes to the output file
    with open(output_dir, "w") as output_file:
        for gene in shared_genes:
            output_file.write(f"{gene}\n")
    shared_count = len(shared_genes)
    print(f"{shared_count} gene names were successfully output into {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="Find shared genes among species assembled gene files")
    parser.add_argument("input_dir", help="Input directory containing species assembled gene files")
    parser.add_argument("species_file", help="Text file containing species names (one per line)")
    parser.add_argument("output_dir", help="Output directory & the output file name for shared genes file")

    args = parser.parse_args()

    # Call the function to find shared genes
    find_shared_genes(args.input_dir, args.species_file, args.output_dir)

if __name__ == "__main__":
    main()
