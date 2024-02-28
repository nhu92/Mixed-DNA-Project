import argparse
from Bio import Phylo
import pandas as pd
import numpy as np

def find_clade_and_move(tree, taxa_name):
    # Find the matching leaf and the smallest clade containing this leaf
    target_leaf = None
    for leaf in tree.get_terminals():  # Iterate through all terminal nodes
        if leaf.name == taxa_name:
            target_leaf = leaf
            break  # Stop when the target leaf is found
    
    # Initialize the list for recording taxa names
    recorded_taxa = []
    
    if target_leaf:
        # Get the path from the root to this leaf, which includes all its ancestors
        path_to_leaf = tree.get_path(target_leaf)
        
        # Iterate from the leaf up to the root
        for clade in reversed(path_to_leaf):
            # Check sister groups of the current clade
            parent = clade
            if parent.clades:  # Ensure this clade has child nodes (i.e., is not a leaf itself)
                for sister in parent.clades:
                    if sister != clade and any("NODE" not in leaf.name for leaf in sister.get_terminals()):
                        # Record names from the sister group if they don't contain "NODE"
                        recorded_taxa.extend([leaf.name for leaf in sister.get_terminals() if "NODE" not in leaf.name])
                        break  # Stop checking other sisters once a valid one is found
                else:
                    # If no valid sister group is found, move up to the upper clade
                    continue
                
                # If the clade has a support value above 0.7, stop the traversal
                if parent.confidence is not None and parent.confidence > 0.7:
                    break
    
    # Return the list of recorded taxa names, excluding any that contain "NODE"
    return recorded_taxa

def calculate_genetic_distance(tree_file):
    # Read the tree file
    tree = tree_file

    # Get the list of terminal node names (i.e., taxa)
    taxa = [leaf.name for leaf in tree.get_terminals()]

    # Initialize a matrix to store the distances
    distance_matrix = np.zeros((len(taxa), len(taxa)))

    # Calculate pairwise distances
    for i in range(len(taxa)):
        for j in range(i+1, len(taxa)):
            distance = tree.distance(taxa[i], taxa[j])
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance

    return taxa, distance_matrix


def genetic_distance_matrix(tree_file, node_output_file, output_file):
    # Load the tree
    tree = Phylo.read(tree_file, 'newick')

    # Reroot logic, assuming this part remains unchanged
    reroot_taxa = ["Amborella", "Nymphaea", "Austrobaileya"]
    for taxa in reroot_taxa:
        for clade in tree.find_clades():
            if clade.name and taxa in clade.name:  # Ensure clade.name is not None
                tree.root_with_outgroup(clade)
                break
        else:
            continue
        break
    else:
        tree.root_at_midpoint()

    # Initialize a list to store records for tips that contain "NODE"
    node_records = []

    # Screen all the tip names (terminal nodes)
    for tip in tree.get_terminals():  # Iterate through terminal nodes only
        if tip.name and "NODE" in tip.name:  # Check if tip has a name and contains "NODE"
            recorded_taxa = find_clade_and_move(tree, tip.name)
            if recorded_taxa:  # Only add if there are recorded taxa
                node_records.append(f'{tip.name}: {"; ".join(recorded_taxa)}')

    # Writing the NODE information to a file
    with open(node_output_file, 'w') as file:
        for record in node_records:
            file.write(f'{record}\n')

    # Proceed with the distance matrix calculation (same as before)
    clades, distance_matrix = calculate_genetic_distance(tree)
    df = pd.DataFrame(distance_matrix)
    df.columns = clades
    df.index = clades
    df.to_csv(output_file)

def main():
    parser = argparse.ArgumentParser(description='Calculate genetic distances and find NODE related taxa in a Newick tree file.')
    parser.add_argument('-t', '--tree', required=True, help='Path to the Newick tree file.')
    parser.add_argument('-n', '--node_output', required=True, help='Path to the output text file for NODE related taxa.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output CSV file for the genetic distance matrix.')
    args = parser.parse_args()

    genetic_distance_matrix(args.tree, args.node_output, args.output)

if __name__ == "__main__":
    main()
