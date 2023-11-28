import argparse
import pandas as pd
import numpy as np
from ete3 import Tree
import matplotlib.pyplot as plt
from skbio.stats.ordination import pcoa

def perform_pcoa(distance_matrix):
    pcoa_results = pcoa(distance_matrix)
    transformed_data = pcoa_results.samples[['PC2', 'PC3']]
    return transformed_data

def plot_pcoa(transformed_data, output_svg):
    plt.figure(figsize=(10, 8))
    plt.scatter(transformed_data['PC2'], transformed_data['PC3'])
    plt.xlabel('PC2')
    plt.ylabel('PC3')
    plt.title('PCoA of Genetic Distances')
    plt.savefig(output_svg, format='svg')

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def sigmoid_trans(data_transformed):
    # Max and min values for scaling, considering non-NaN values only
    max_value = data_transformed.max().max()
    min_value = data_transformed[data_transformed != 0].min().min()

    # Scaling the data for sigmoid transformation
    scaled_data = (data_transformed - min_value) / (max_value - min_value)

    # Apply sigmoid transformation
    sigmoid_transformed = sigmoid(scaled_data * 10 - 5)  # Scale and shift for effective sigmoid transformation

    # Ensuring diagonal elements are zero
    np.fill_diagonal(sigmoid_transformed.values, 0)

    return sigmoid_transformed

def main():
    parser = argparse.ArgumentParser(description='PCoA on Phylogenetic Tree')
    parser.add_argument('input_matrix', type=str, help='Input distance matrix file')
    parser.add_argument('output_csv', type=str, help='Output CSV file for PCoA data')
    parser.add_argument('output_svg', type=str, help='Output SVG file for plot')
    args = parser.parse_args()

    distance_matrix = pd.read_csv(args.input_matrix, header=0, index_col=0)
    siged_matrix = sigmoid_trans(distance_matrix)
    transformed_data = perform_pcoa(siged_matrix)
    transformed_data.index = distance_matrix.index  # Add species names as index
    transformed_data.to_csv(args.output_csv, index=True)
    plot_pcoa(transformed_data, args.output_svg)

if __name__ == "__main__":
    main()

