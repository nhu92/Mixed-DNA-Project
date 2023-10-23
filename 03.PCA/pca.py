import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import argparse
import matplotlib.colors as mcolors

def pca_analysis(matrix_file, table_file, output_file):
    # Load the matrix and table from the input files
    matrix = pd.read_csv(matrix_file, index_col=0)
    table = pd.read_csv(table_file)
    matrix.columns = matrix.columns.str.replace(r'\_mean', '', regex=True)
    matrix = matrix.T

    # Normalize the matrix
    scaler = StandardScaler()
    normalized_matrix = scaler.fit_transform(matrix)
    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(normalized_matrix)

    # Create a DataFrame for the PCA result
    pca_df = pd.DataFrame(data = pca_result, columns = ['PCA1', 'PCA2'], index=matrix.index)

    # Create a color map and taxa map from the table
    color_map = {}
    taxa_map = {}
    for index, row in table.iterrows():
        color_map[row['item']] = row['color']
        taxa_map[row['item']] = row['taxa']

    # Assign colors and taxa to each dot in the PCA plot
    colors = []
    taxa = []
    for item in matrix.index:
        color_assigned = False
        for table_item, color in color_map.items():
            if item.startswith(table_item):
                rgba_color = mcolors.to_rgba(color, alpha=0.6)  # Convert the color name to RGBA and make it 60% transparent
                colors.append(rgba_color)
                taxa.append(taxa_map[table_item])
                color_assigned = True
                break
        if not color_assigned:
            colors.append('red')
            taxa.append('Undetermined Species')

    # Create the PCA plot
    fig, ax = plt.subplots(figsize=(8,8))
    scatter = ax.scatter(pca_df['PCA1'], pca_df['PCA2'], c=colors)
    
    # Show the explained variance ratio on the axes
    ax.set_xlabel(f'Principal Component 1 ({pca.explained_variance_ratio_[0]*100:.2f}%)')
    ax.set_ylabel(f'Principal Component 2 ({pca.explained_variance_ratio_[1]*100:.2f}%)')
    
    # Create a legend for the plot using the taxa information
    legend_labels = dict(zip(taxa, colors))
    
    # Order the legend labels alphabetically
    sorted_legend_labels = dict(sorted(legend_labels.items()))
    
    handles = [plt.Line2D([], [], marker='o', color=color, linestyle='', markersize=10) for color in sorted_legend_labels.values()]
    
    ax.legend(handles, sorted_legend_labels.keys(), title='Taxa', loc='upper right')
    
    plt.title('2 Component PCA')
    
    # Save the plot as an SVG file
    plt.savefig(output_file, format='svg')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform PCA analysis.')
    parser.add_argument('matrix', help='Path to the input matrix file.')
    parser.add_argument('table', help='Path to the input table file.')
    parser.add_argument('output', help='Path to the output SVG file.')
    
    args = parser.parse_args()
    
    pca_analysis(args.matrix, args.table, args.output)
