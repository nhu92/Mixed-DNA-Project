import argparse
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
import re

def calculate_euclidean_distance(coordinates):
    distance_matrix = pdist(coordinates, 'euclidean')
    distance_matrix = squareform(distance_matrix)
    np.fill_diagonal(distance_matrix, np.inf)
    return distance_matrix

def find_best_matches(distance_matrix, nodes, row_names, num_matches=3):
    best_matches = []

    for node in nodes:
        matches = []
        for row_idx, row in enumerate(distance_matrix[node]):
            if row_idx not in nodes:
                matches.append((row_names[row_idx], row))
        matches = sorted(matches, key=lambda x: x[1])[:num_matches]
        best_matches.append((row_names[node], [match[0] for match in matches]))

    return best_matches

def remove_pattern(text):
    return re.sub(r'_\d{4}_mean', '', text)

def main():
    parser = argparse.ArgumentParser(description='Find the best matches for 3D coordinates.')
    parser.add_argument('input_file', help='Input CSV file with 3D coordinates.')
    parser.add_argument('output_file', help='Output CSV file with the best matches.')
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

    # Find the best matches for nodes (excluding the "NODE" row)
    best_matches = find_best_matches(distance_matrix, nodes, row_names)

    # Create a new DataFrame with the best matches
    output_data = []
    for node, matches in best_matches:
        for match in matches:
            output_data.append((remove_pattern(node), remove_pattern(match)))
    
    output_df = pd.DataFrame(output_data, columns=['Node', 'BestMatches'])
    output_df.to_csv(args.output_file, index=False)

if __name__ == '__main__':
    main()
