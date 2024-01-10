import argparse
from Bio import Phylo
import pandas as pd
import numpy as np

def calculate_genetic_distance(tree_file):
    # Read the tree file and assume you have a Bio.Phylo tree object

    # Get the list of terminal node names (i.e., taxa)
    taxa = [leaf.name for leaf in tree_file.get_terminals()]

    # Initialize a matrix to store the distances
    distance_matrix = np.zeros((len(taxa), len(taxa)))

    # Calculate pairwise distances
    for i in range(len(taxa)):
        for j in range(i+1, len(taxa)):
            try:
                # Get the common ancestor of taxa[i] and taxa[j]
                common_ancestor = tree_file.common_ancestor(taxa[i], taxa[j])
                # Calculate the patristic distance (sum of branch lengths in the path)
                distance = common_ancestor.get_distance(tree_file)
                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance
            except:
                # Handle cases where the taxa pair is not in the tree
                pass

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
