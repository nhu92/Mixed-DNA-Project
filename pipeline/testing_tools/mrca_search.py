import csv
from Bio import Phylo
import argparse

def reroot_tree(tree):
    # Reroot the tree by 'amborella', if not exist then 'water lily', if both not exist, root by midpoint.
    for root_name in ['amborella', 'water_lily']:
        root_clades = list(tree.find_clades(name=root_name))
        if root_clades:
            tree.root_with_outgroup(root_clades[0])
            return tree
    tree.root_at_midpoint()
    return tree

def analyze_unknown_species(tree):
    results = {}
    # Find unknown species (any species with 'NODE' in its name).
    unknown_species = [clade for clade in tree.get_terminals() if 'NODE' in clade.name]
    
    for unknown in unknown_species:
        path_to_root = tree.get_path(unknown)
        for clade in reversed(path_to_root):  # Start from the node closest to the unknown and move up
            if clade.confidence is not None and clade.confidence > 0.75:
                # Get all species or clades within this clade, excluding unknowns
                clade_group = [term.name for term in clade.get_terminals() if 'NODE' not in term.name]
                results[unknown.name] = clade_group
                break  # Stop at the first clade with support > 0.75
    return results

def write_to_csv(output_file, analysis_results):
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['Unknown_Species', 'Related_Clad']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for unknown, clade_group in analysis_results.items():
            writer.writerow({'Unknown_Species': unknown, 'Related_Clade': ', '.join(clade_group)})

def main(args):
    tree = Phylo.read(args.input, 'newick')
    tree = reroot_tree(tree)
    analysis_results = analyze_unknown_species(tree)
    
    # Write the analysis results to a CSV file
    write_to_csv(args.output, analysis_results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze Unknown Species in Phylogenetic Tree.')
    parser.add_argument('-i', '--input', type=str, help='Input file path for the Newick tree', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output CSV file path', required=True)
    args = parser.parse_args()
    main(args)
