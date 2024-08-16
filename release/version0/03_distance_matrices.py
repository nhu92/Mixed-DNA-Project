import subprocess
import os
import argparse
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor
from Bio import Phylo
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

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

def group_and_sum(input_file, output_file):
    data = pd.read_csv(input_file)
    filtered_data = data[~data['row_name'].str.contains("NODE")]
    grouped_sum = filtered_data.groupby('row_name')['total_value'].sum().reset_index()
    grouped_sum.to_csv(output_file, index=False)

def log_status(log_file, message):
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(log_file, 'a') as log:
        log.write(f"[{timestamp}] {message}\n")
        log.flush()

def run_command(command, step_name, log_file):
    try:
        subprocess.run(command, shell=True, check=True)
        log_status(log_file, f"{step_name}: SUCCESS")
    except subprocess.CalledProcessError:
        log_status(log_file, f"{step_name}: FAILURE")
        print(f"Error: {step_name} failed. Check {log_file} for details.")
        exit(1)

def process_gene(gene_name_shorter, tree_dir, log_file):
    try:
        tree_list_file = "./loop.treelist.txt"
        run_command(f'ls "{tree_dir}/{gene_name_shorter}"*"tre" > {tree_list_file}', f"List trees for {gene_name_shorter}", log_file)

        with open(tree_list_file, 'r') as tree_files:
            for i, filename in enumerate(tree_files, start=1):
                filename = filename.strip()

                node_output_file = f"{tree_dir}/{gene_name_shorter}.{i}.list.txt"
                output_file = f"{tree_dir}/{gene_name_shorter}.{i}.matrix"
                
                genetic_distance_matrix(filename, node_output_file, output_file)
                log_status(log_file, f"Generated matrix for {gene_name_shorter} tree {i}")

                copy_cmd = f'cp "{output_file}" "{tree_dir}/{gene_name_shorter}.{i}.cleaned.csv"'
                run_command(copy_cmd, f"Copy matrix to cleaned CSV for {gene_name_shorter} tree {i}", log_file)

        os.remove(tree_list_file)
        log_status(log_file, "Removed temporary file loop.treelist.txt")
    except Exception as e:
        log_status(log_file, f"Failed processing {gene_name_shorter}: {e}")
        print(f"Failed processing {gene_name_shorter}: {e}")

def main(threads, proj_name, gene_list, log_file, threshold, use_flag):
    tree_dir = "04_all_trees"
    os.makedirs(tree_dir, exist_ok=True)
    log_status(log_file, "Created directory ./04_all_trees")

    run_command(f'cp 03_phylo_results/*.tre {tree_dir}', "Copy tree files to ./04_all_trees", log_file)

    with open(gene_list, 'r') as genes:
        gene_names = [gene.strip() for gene in genes]

    with ThreadPoolExecutor(max_workers=int(threads)) as executor:
        for gene_name_shorter in gene_names:
            executor.submit(process_gene, gene_name_shorter, tree_dir, log_file)

    final_output = process_matrices(tree_dir, proj_name, threshold, use_flag)
    output_file = f'./{proj_name}.summary_dist.csv'
    final_output.to_csv(output_file, index=False)
    log_status(log_file, f"Processed matrices saved to {output_file}")

    # Run group and sum directly
    cumulative_output_file = f'./{proj_name}.cumulative_dist.csv'
    group_and_sum(output_file, cumulative_output_file)
    log_status(log_file, f"Generated cumulative distance file: {cumulative_output_file}")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline for distance calculation and tree processing.")
    parser.add_argument("-t", "--threads", required=True, help="Number of threads to use for parallel processing")
    parser.add_argument("-p", "--proj_name", required=True, help="Project name")
    parser.add_argument("-g", "--gene_list", type=str, default="gene_list.txt", help="Path to the gene list file")
    parser.add_argument("--threshold", type=float, default=1.96, help="Threshold for value adjustment")
    parser.add_argument("--use_flag", action="store_true", help="Use flag to assign the smallest number as 0 and others as 999")

    args = parser.parse_args()

    log_file = f"{args.proj_name}_03_distance_calc.log"

    log_status(log_file, "Pipeline started with the following parameters:")
    log_status(log_file, f"  Threads: {args.threads}")
    log_status(log_file, f"  Project Name: {args.proj_name}")
    log_status(log_file, f"  Gene List: {args.gene_list}")
    log_status(log_file, f"  Threshold: {args.threshold}")
    log_status(log_file, f"  Use Flag: {args.use_flag}")

    main(args.threads, args.proj_name, args.gene_list, log_file, args.threshold, args.use_flag)

    log_status(log_file, "Pipeline completed successfully.")
    print(f"Pipeline completed. Check {log_file} for the status of each step.")