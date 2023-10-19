import argparse
from Bio import Phylo
import pandas as pd

def calculate_relative_distance(tree_file, node_str, output_file):
    # Load the tree
    tree = Phylo.read(tree_file, 'newick')

    # Find the nodes
    nodes = [clade for clade in tree.find_clades() if clade.name and node_str in clade.name]

    if not nodes:
        print(f"No nodes found with string {node_str} in their names.")
        return

    # Calculate the distance from each node to all other nodes
    distances = {}
    for node in nodes:
        for clade in tree.find_clades():
            if clade.name not in distances:
                distances[clade.name] = {}
            distances[clade.name][node.name] = tree.distance(node, clade)

    # Output to a table
    df = pd.DataFrame(distances).T
    df.to_csv(output_file)

def main():
    parser = argparse.ArgumentParser(description='Calculate relative distance between nodes in a Newick tree file.')
    parser.add_argument('-t', '--tree', required=True, help='Path to the Newick tree file.')
    parser.add_argument('-n', '--node', required=True, help='String to search for in node names.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output CSV file.')
    args = parser.parse_args()

    calculate_relative_distance(args.tree, args.node, args.output)

if __name__ == "__main__":
    main()
