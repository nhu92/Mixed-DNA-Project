import argparse
from ete3 import Tree

def collapse_repeated_tips(tree):
    leaf_names = [leaf.name for leaf in tree.iter_leaves()]
    for name in leaf_names:
        leaves = [leaf for leaf in tree.iter_leaves() if leaf.name == name]
        if len(leaves) > 1:
            common_ancestor = tree.get_common_ancestor(leaves)
            for leaf in leaves:
                leaf.detach()
            common_ancestor.add_child(name=name)

def are_same_topology(tree1, tree2):
    collapse_repeated_tips(tree1)
    collapse_repeated_tips(tree2)
    tree1_sorted = tree1.get_topology_id(attr="name")
    tree2_sorted = tree2.get_topology_id(attr="name")
    return tree1_sorted == tree2_sorted

def read_tree(filename):
    tree = Tree(filename, format=1)
    for node in tree.traverse():
        node.dist = 0.0
        node.name = ""
    return tree

def main():
    parser = argparse.ArgumentParser(description='Compare the topology of two phylogenetic trees.')
    parser.add_argument('tree1', type=str, help='The filename of the first tree.')
    parser.add_argument('tree2', type=str, help='The filename of the second tree.')
    args = parser.parse_args()

    # Read trees from files
    tree1 = read_tree(args.tree1)
    tree2 = read_tree(args.tree2)

    # Compare the topologies of the trees
    print(are_same_topology(tree1, tree2))

if __name__ == "__main__":
    main()
