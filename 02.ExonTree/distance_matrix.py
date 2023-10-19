import argparse
from Bio import Phylo
import pandas as pd
import numpy as np

def calculate_genetic_distance(tree_file):
    # Read the tree file
    tree = tree_file

    # Get the list of terminal node names (i.e., taxa)
    taxa = [leaf.name for leaf in tree.get_terminals()]

    # Initialize a matrix to store the distances
    distance_matrix = np.zeros((len(taxa), len(taxa)))

    # Calculate pairwise distances
    for i in range(len(taxa)):
        for j in range(i+1, len(taxa)):
            distance = tree.distance(taxa[i], taxa[j])
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance

    return taxa, distance_matrix


def genetic_distance_matrix(tree_file, output_file):
    # Load the tree
    tree = Phylo.read(tree_file, 'newick')

    clades, distance_matrix = calculate_genetic_distance(tree)

    df = pd.DataFrame(distance_matrix)
    df.columns = clades
    df.index = clades
    # Output to a table
    df.to_csv(output_file)

def main():
    parser = argparse.ArgumentParser(description='Calculate relative distance between nodes in a Newick tree file.')
    parser.add_argument('-t', '--tree', required=True, help='Path to the Newick tree file.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output CSV file.')
    args = parser.parse_args()

    genetic_distance_matrix(args.tree, args.output)

if __name__ == "__main__":
    main()
