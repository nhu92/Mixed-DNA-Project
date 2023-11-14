import argparse
import pandas as pd
import numpy as np
from ete3 import Tree
import matplotlib.pyplot as plt
from skbio.stats.ordination import pcoa

# PC1-PC3
def perform_pcoa(distance_matrix):
    pcoa_results = pcoa(distance_matrix)
    transformed_data = pcoa_results.samples[['PC1', 'PC3']]
    return transformed_data

def plot_pcoa(transformed_data, output_svg):
    plt.figure(figsize=(10, 8))
    plt.scatter(transformed_data['PC1'], transformed_data['PC3'])
    plt.xlabel('PC1')
    plt.ylabel('PC3')
    plt.title('PCoA of Genetic Distances')
    plt.savefig(output_svg, format='svg')

def main():
    parser = argparse.ArgumentParser(description='PCoA on Phylogenetic Tree')
    parser.add_argument('input_matrix', type=str, help='Input distance matrix file')
    parser.add_argument('output_csv', type=str, help='Output CSV file for PCoA data')
    parser.add_argument('output_svg', type=str, help='Output SVG file for plot')
    args = parser.parse_args()

    distance_matrix = pd.read_csv(args.input_matrix, header=0, index_col=0)
    transformed_data = perform_pcoa(distance_matrix)
    transformed_data.index = distance_matrix.index  # Add species names as index
    transformed_data.to_csv(args.output_csv, index=True)
    plot_pcoa(transformed_data, args.output_svg)

if __name__ == "__main__":
    main()
