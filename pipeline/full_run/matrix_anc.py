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
        for j in range(i + 1, len(taxa)):
            # Find the last common ancestor
            lca = tree.common_ancestor({'name': taxa[i]}, {'name': taxa[j]})

            # Function to get distance to parent or LCA if directly connected
            def get_distance_to_parent_or_lca(taxon):
                path = tree.get_path(taxon)
                if len(path) > 1:
                    return tree.distance(lca, path[-2])  # Distance to parent
                else:
                    return 0  # Taxon is directly connected to the LCA

            # Calculate the distances
            distance_to_parent_i = get_distance_to_parent_or_lca(taxa[i])
            distance_to_parent_j = get_distance_to_parent_or_lca(taxa[j])

            # Calculate the genetic distance
            genetic_distance = distance_to_parent_i + distance_to_parent_j
            distance_matrix[i, j] = genetic_distance
            distance_matrix[j, i] = genetic_distance

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
