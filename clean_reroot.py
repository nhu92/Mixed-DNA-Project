import argparse
from ete3 import Tree

def trim_leaves(tree, start_str):
    for leaf in tree.iter_leaves():
        if leaf.name.startswith(start_str):
            leaf.delete()

def reroot_tree(tree, names):
    for name in names:
        leaves = tree.get_leaves_by_name(name)
        if leaves:
            tree.set_outgroup(leaves[0])
            return
    tree.set_outgroup(tree.get_midpoint_outgroup())

def main():
    parser = argparse.ArgumentParser(description='Process a tree.')
    parser.add_argument('--tree', type=str, required=True, help='The input tree file.')
    parser.add_argument('--start_str', type=str, required=True, help='The input string.')
    parser.add_argument('--output', type=str, required=True, help='The output file.')
    args = parser.parse_args()

    with open(args.tree, 'r') as f:
        trees = f.read().split(';\n')

    output_trees = []
    names = [f'amborella_{i}' for i in range(10000)] + [f'water_lily_{i}' for i in range(10000)]
    
    for t in trees:
        if t.strip() == '':
            continue
        tree = Tree(t + ';')
        trim_leaves(tree, args.start_str)
        reroot_tree(tree, names)
        output_trees.append(tree.write(format=1))

    with open(args.output, 'w') as f:
        f.write(';\n'.join(output_trees) + ';')

if __name__ == "__main__":
    main()
