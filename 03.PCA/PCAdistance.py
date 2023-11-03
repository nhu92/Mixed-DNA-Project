import argparse
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform

def calculate_euclidean_distance(coordinates):
    distance_matrix = pdist(coordinates, 'euclidean')
    distance_matrix = squareform(distance_matrix)
    np.fill_diagonal(distance_matrix, np.inf)
    return distance_matrix

def find_best_matches(distance_matrix, nodes, row_names, k=3):
    best_matches = []

    for node in nodes:
        distances = distance_matrix[node]
        match_indices = np.argpartition(distances, k)[:k]
        best_matches.append([row_names[node]] + [row_names[i] for i in match_indices])

    return best_matches

def main():
    parser = argparse.ArgumentParser(description='Find the best 3 distance matches for 3D coordinates.')
    parser.add_argument('input_file', help='Input CSV file with 3D coordinates.')
    parser.add_argument('output_file', help='Output CSV file with best 3 distance matches.')
    args = parser.parse_args()

    # Read the input CSV file
    df = pd.read_csv(args.input_file)

    # Extract coordinates and row names
    coordinates = df.iloc[:, 1:].values
    row_names = df.iloc[:, 0].tolist()

    # Define the nodes that start with "NODE"
    nodes = [i for i, name in enumerate(row_names) if name.startswith("NODE")]

    # Calculate pairwise Euclidean distances
    distance_matrix = calculate_euclidean_distance(coordinates)

    # Find the best 3 distance matches for nodes
    best_matches = find_best_matches(distance_matrix, nodes, row_names, k=3)

    # Create a new DataFrame with the best 3 distance matches
    output_df = pd.DataFrame(best_matches, columns=['Node'] + [f'Match_{i}' for i in range(1, 4)])
    output_df.to_csv(args.output_file, index=False)

if __name__ == '__main__':
    main()
