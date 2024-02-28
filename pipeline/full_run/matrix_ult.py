import argparse
from Bio import Phylo
import pandas as pd
import numpy as np

def find_clade_and_move(tree, taxa_name):
    # Find the clade by taxa name
    target_clade = None
    for clade in tree.find_clades():
        if clade.name == taxa_name:
            target_clade = clade
            break
    
    recorded_taxa = []
    while target_clade:
        # Check sister groups
        parent = tree.get_path(target_clade)[-2] if len(tree.get_path(target_clade)) > 1 else None
        if parent:
            for sister in parent.clades:
                if sister != target_clade and any("NODE" not in leaf.name for leaf in sister.get_terminals()):
                    recorded_taxa.append(target_clade.name)
                    break
            else:
                # If no valid sister, move to upper clade
                target_clade = parent
                continue

            # Move to upper clade and check support value
            target_clade = parent
            if parent.confidence is not None and parent.confidence > 0.7:
                return [taxa for taxa in recorded_taxa if "NODE" not in taxa]
        break  # Exit while loop if there's no parent or if other conditions are not met

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

    # Reroot logic
    reroot_taxa = ["Amborella", "Nymphaea", "Austrobaileya"]
    for taxa in reroot_taxa:
        for clade in tree.find_clades():
            if clade.name and taxa in clade.name:  # Check if clade.name is not None and contains taxa
                tree.root_with_outgroup(clade)
                break
        else:
            continue  # Continue if the taxa was not found in any clade
        break  # Break if the tree was rerooted successfully
    else:
        # If none of the specified taxa were found, reroot at the midpoint
        tree.root_at_midpoint()

    # Screen all the taxa names containing "NODE"
    node_related_taxa = {}
    for clade in tree.find_clades():
        if "NODE" in clade.name:
            associated_taxa = find_clade_and_move(tree, clade.name)
            if associated_taxa:
                node_related_taxa[clade.name] = associated_taxa

    # Save the NODE related taxa into a file
    with open(node_output_file, 'w') as f:
        for node, taxa_list in node_related_taxa.items():
            f.write(f"{node}: {', '.join(taxa_list)}\n")

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
