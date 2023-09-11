import os

# Function to get the shared genes among all files
def find_shared_genes(input_dir, species_names, output_dir):
    # Initialize a dictionary to store genes for each species
    species_genes = {species: set() for species in species_names}

    # List files in the input directory
    files = os.listdir(input_dir)

    # Iterate through files and extract genes
    for file_name in files:
        if file_name.endswith("_assembled_genes.txt"):
            species_name = file_name.split("_")[0]
            file_path = os.path.join(input_dir, file_name)

            with open(file_path, "r") as file:
                for line in file:
                    gene_name = line.strip()
                    species_genes[species_name].add(gene_name)

    # Find the shared genes among all species
    shared_genes = set.intersection(*species_genes.values())

    # Write the shared genes to the output file
    output_path = os.path.join(output_dir, "shared_genes.txt")
    with open(output_path, "w") as output_file:
        for gene in shared_genes:
            output_file.write(f"{gene}\n")

# List of species names (modify as needed)
species_names = ["species1", "species2", "species3"]

# Input directory containing files
input_directory = "input_directory_path"  # Modify this to your input directory path

# Output directory (create if it doesn't exist)
output_directory = "output_directory_path"  # Modify this to your output directory path
os.makedirs(output_directory, exist_ok=True)

# Call the function to find shared genes
find_shared_genes(input_directory, species_names, output_directory)
