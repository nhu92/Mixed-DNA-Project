import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from ete3 import Tree

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Generate a heatmap with a phylogenetic tree.')
parser.add_argument('data_file', type=str, help='Path to the data file (CSV format).')
parser.add_argument('tree_file', type=str, help='Path to the phylogenetic tree file (Newick format).')
parser.add_argument('output_file', type=str, help='Path to the output file (SVG format).')
args = parser.parse_args()

# Load your data
data = pd.read_csv(args.data_file, index_col=0)

# Separate the mean and standard deviation columns
mean_columns = [col for col in data.columns if 'mean' in col]
sd_columns = [col for col in data.columns if 'sd' in col]
mean_data = data[mean_columns]
sd_data = data[sd_columns]

# Create a dendrogram based on the phylogenetic tree
tree = Tree(args.tree_file)
leaves = tree.get_leaf_names()
linkage_matrix = linkage(pdist(mean_data.loc[leaves].values), method='average')

# Create a figure with two subplots: one for the dendrogram (tree) and one for the heatmap
fig = plt.figure(figsize=(10, 8))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3]) 

# Create a dendrogram (tree)
ax0 = plt.subplot(gs[0])
dendrogram(linkage_matrix, orientation='left', labels=leaves, ax=ax0)
ax0.set_xticks([])
ax0.set_yticks([])
ax0.axis('off')  # Remove the frame

# Create a heatmap
ax1 = plt.subplot(gs[1])
sns.heatmap(mean_data.loc[leaves], cmap='RdBu_r', annot=sd_data.loc[leaves], fmt=".2f", ax=ax1)

# Rename the x-axis labels
new_labels = [label.replace('_mean', '') for label in mean_data.columns]
ax1.set_xticklabels(new_labels)

plt.savefig(args.output_file)
