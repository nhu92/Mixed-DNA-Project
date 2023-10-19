import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import argparse

def pca_analysis(matrix_file, table_file, output_file):
    # Load the matrix and table from the input files
    matrix = pd.read_csv(matrix_file, index_col=0)
    table = pd.read_csv(table_file)

    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(matrix)

    # Create a DataFrame for the PCA result
    pca_df = pd.DataFrame(data = pca_result, columns = ['PCA1', 'PCA2'], index=matrix.index)

    # Create a color map from the table
    color_map = {}
    for index, row in table.iterrows():
        color_map[row['item']] = row['color']

    # Assign colors to each dot in the PCA plot
    colors = []
    for item in matrix.index:
        color_assigned = False
        for table_item, color in color_map.items():
            if item.startswith(table_item):
                colors.append(color)
                color_assigned = True
                break
        if not color_assigned:
            colors.append('red')

    # Create the PCA plot
    plt.figure(figsize=(8,8))
    plt.scatter(pca_df['PCA1'], pca_df['PCA2'], c=colors)
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('2 Component PCA')
    
    # Save the plot as an SVG file
    plt.savefig(output_file, format='svg')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform PCA analysis.')
    parser.add_argument('--matrix', required=True, help='Path to the input matrix file.')
    parser.add_argument('--table', required=True, help='Path to the input table file.')
    parser.add_argument('--output', required=True, help='Path to the output SVG file.')
    
    args = parser.parse_args()
    
    pca_analysis(args.matrix, args.table, args.output)
