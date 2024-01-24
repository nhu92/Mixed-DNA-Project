import argparse
from Bio import Phylo
import numpy as np
import pandas as pd

def pairwise_distances(trees, species):
    """Calculate pairwise distances for all species across all trees."""
    # Initialize a dictionary to hold the sum of distances and counts
    distances_sum = {(s1, s2): 0.0 for s1 in species for s2 in species if s1 != s2}
    counts = {(s1, s2): 0 for s1 in species for s2 in species if s1 != s2}

    for tree in trees:
        for s1 in species:
            for s2 in species:
                if s1 != s2:
                    try:
                        # Get the common ancestor of s1 and s2
                        common_ancestor = tree.common_ancestor(s1, s2)
                        # Calculate the patristic distance (sum of branch lengths in the path)
                        distance = common_ancestor.get_distance(tree)
                        distances_sum[(s1, s2)] += distance
                        counts[(s1, s2)] += 1
                    except:
                        # Handle cases where the species pair is not in the tree
                        pass

    # Average the distances
    for pair in distances_sum:
        if counts[pair] > 0:
            distances_sum[pair] /= counts[pair]

    return distances_sum


def create_distance_matrix(distances, species):
    """Create a distance matrix from the averaged distances."""
    matrix = np.zeros((len(species), len(species)))

    for i, s1 in enumerate(species):
        for j, s2 in enumerate(species):
            if i != j:
                matrix[i, j] = distances.get((s1, s2), 0.0)
                matrix[j, i] = matrix[i, j]  # Symmetric matrix

    return pd.DataFrame(matrix, index=species, columns=species)

def main():
    parser = argparse.ArgumentParser(description='Create Averaged Pairwise Distance Matrix from Gene Trees')
    parser.add_argument('input_file', type=str, help='Input Newick tree file')
    parser.add_argument('output_file', type=str, help='Output file for the distance matrix')
    args = parser.parse_args()

    trees = list(Phylo.parse(args.input_file, 'newick'))
    species = set()
    for tree in trees:
        species.update([term.name for term in tree.get_terminals()])
    
    # Convert species set to a list for consistent ordering
    species_list = sorted(list(species))

    distances = pairwise_distances(trees, species_list)
    distance_matrix = create_distance_matrix(distances, species_list)

    distance_matrix.to_csv(args.output_file, sep=',')

if __name__ == "__main__":
    main()


