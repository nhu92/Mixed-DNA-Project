import argparse
import re
from ete3 import Tree
from Bio import Phylo
import os
import sys

def check_file_exists(filename):
    if not os.path.isfile(filename):
        print(f"Error: The file '{filename}' does not exist.")
        sys.exit()


def add_leaf_next_to_matching_tip(tree, leaf_name):
    matching_tip = tree.search_nodes(name=leaf_name)

    if matching_tip:
        # Get the matching tip
        matching_tip = matching_tip[0]
        # Get the parent of the matching tip
        parent = matching_tip.up
        # Create a new clade and add it as a sibling to the matching tip
        new_clade = parent.add_child(name=f"Clade_{leaf_name}")
        # Move the matching tip and create a new leaf inside the new clade
        matching_tip.detach()
        new_clade.add_child(matching_tip)
        new_leaf = new_clade.add_child(name=leaf_name)
    else:
        print(f"No matching tip found in tree2 for {leaf_name}")

def expand_leaves(tree1, tree2):
    # Get the list of leaf names in both trees
    leaf_names1 = sorted([leaf.name for leaf in tree1.iter_leaves()])
    leaf_names2 = sorted([leaf.name for leaf in tree2.iter_leaves()])

    # Create a list of leaves to add to tree2
    list1 = [item for item in leaf_names2 if leaf_names1.count(item) > 1]
    # Traverse list1 and tree2 to add leaves next to matching tips
    for leaf_name in list1:
        add_leaf_next_to_matching_tip(tree2, leaf_name)

def are_same_topology(tree1, tree2):
    # Unroot the trees
    tree1.unroot()
    tree2.unroot()

    # Expand leaves in tree2
    expand_leaves(tree1, tree2)

    # Sort the leaf node names and compare
    tree1_sorted = tree1.get_topology_id(attr="name")
    tree2_sorted = tree2.get_topology_id(attr="name")
    print(tree1_sorted)
    print(tree2_sorted)
    return tree1_sorted == tree2_sorted

def remove_patterns_from_newick(newick_str):
    # Remove all instances of '_R_'
    newick_str = re.sub(r'_R_', '', newick_str)

    # Remove all instances of 'numbers_'
    newick_str = re.sub(r'\d+_', '', newick_str)

    # Remove the last '_number' pattern
    newick_str = re.sub(r'_\d+', '', newick_str)

    return newick_str

def read_tree(filename):
    # Read the Newick string from the file
    with open(filename, 'r') as file:
        newick_str = file.read()

    # Remove certain patterns from the Newick string
    newick_str = remove_patterns_from_newick(newick_str)

    # Read the tree from the Newick string
    tree = Tree(newick_str, format=1)

    # Remove all the branch lengths, node labels, etc.
    #for node in tree.traverse():
    #    node.dist = 0.0  # Remove branch length
    #    node.name = ""  # Remove node label

    return tree

def main():
    parser = argparse.ArgumentParser(description='Compare the topology of two phylogenetic trees.')
    parser.add_argument('tree1', type=str, help='The filename of the first tree.')
    parser.add_argument('tree2', type=str, help='The filename of the second tree.')
    args = parser.parse_args()

    # Read trees from files
    tree1 = read_tree(args.tree1)
    tree2 = read_tree(args.tree2)

    check_file_exists(args.tree1)
    # Compare the topologies of the trees
    print(args.tree1, "matches target topology:", are_same_topology(tree1, tree2))

if __name__ == "__main__":
    main()
