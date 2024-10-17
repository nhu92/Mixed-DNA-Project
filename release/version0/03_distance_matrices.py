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
    for leaf in tree.get_terminals():
        if leaf.name == taxa_name:
            target_leaf = leaf
            break
    recorded_taxa = []
    
    if target_leaf:
        path_to_leaf = tree.get_path(target_leaf)
        for clade in reversed(path_to_leaf):
            parent = clade
            if parent.clades:
                for sister in parent.clades:
                    if sister != clade and any("NODE" not in leaf.name for leaf in sister.get_terminals()):
                        recorded_taxa.extend([leaf.name for leaf in sister.get_terminals() if "NODE" not in leaf.name])
                        break
                else:
                    continue
                if parent.confidence is not None and parent.confidence > 0.7:
                    break
    return recorded_taxa

def calculate_genetic_distance(tree_file):
    tree = tree_file
    taxa = [leaf.name for leaf in tree.get_terminals()]
    distance_matrix = np.zeros((len(taxa), len(taxa)))
    for i in range(len(taxa)):
        for j in range(i+1, len(taxa)):
            distance = tree.distance(taxa[i], taxa[j])
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance
    return taxa, distance_matrix

def distance_to_similarity(distance_df):
    similarity_df = distance_df.copy()
    numeric_cols = similarity_df.select_dtypes(include=['float64', 'int']).columns
    similarity_df[numeric_cols] = 1 / (1 + similarity_df[numeric_cols])
    return similarity_df

def clean_up_matrix(df, proj_name, threshold, taxa_file, use_flag=False):
    df = df[~df[df.columns[0]].str.contains(proj_name)]
    cols_to_keep = [df.columns[0]] + [col for col in df.columns[1:] if proj_name in col]
    df = df[cols_to_keep]
    for col in df.columns[1:]:
        if use_flag:
            min_value = df[col].min()
            df[col] = df[col].apply(lambda x: 0 if x == min_value else 999)
        else:
            col_mean = df[col].mean()
            col_sd = df[col].std()
            df[col] = df[col].apply(lambda x: 999 if x > (col_mean - threshold * col_sd) else x)

        if taxa_file:
            species_to_taxa = {}
            with open(taxa_file, 'r') as file:
                for line in file:
                    parts = line.strip().split(':')
                    if len(parts) == 2:
                        species, taxa = parts
                        species_to_taxa[species.strip()] = [tax.strip() for tax in taxa.split(';')]
            for header in df.columns[1:]:
                for species, taxa_list in species_to_taxa.items():
                    if species == header:
                        for index, row in df.iterrows():
                            if row[df.columns[0]] not in taxa_list:
                                df.at[index, header] = 999 
    df.iloc[:, 0] = df.iloc[:, 0].str.replace(r'\d+', '', regex=True).str.rstrip('_')
    return df

def process_matrices(directory, proj_name, threshold, flag):
    all_matrices = []
    for filename in os.listdir(directory):
        if filename.endswith('cleaned.csv'):
            matrix_path = os.path.join(directory, filename)
            matrix = pd.read_csv(matrix_path)
            prefix = filename.split("cleaned.csv")[0]
            list_file = f"{prefix}list.txt"
            list_file_path = os.path.join(directory, list_file)
            matrix = clean_up_matrix(matrix, proj_name, threshold, list_file_path, flag)
            matrix = distance_to_similarity(matrix)
            all_matrices.append(matrix)
    total_matrix = pd.concat(all_matrices, ignore_index=True)
    total_matrix.fillna(0, inplace=True)
    numeric_cols = [col for col in total_matrix.columns if col != 'Unnamed: 0']
    total_matrix['total_value'] = total_matrix[numeric_cols].sum(axis=1)
    final_output = total_matrix[['Unnamed: 0', 'total_value']].rename(columns={'Unnamed: 0': 'row_name'})
    return final_output

def genetic_distance_matrix(tree_file, node_output_file, output_file):
    tree = Phylo.read(tree_file, 'newick')
    reroot_taxa = ["Amborella", "Nymphaea", "Austrobaileya"]
    for taxa in reroot_taxa:
        for clade in tree.find_clades():
            if clade.name and taxa in clade.name:
                tree.root_with_outgroup(clade)
                break
        else:
            continue
        break
    else:
        tree.root_at_midpoint()
    node_records = []
    for tip in tree.get_terminals():
        if tip.name and "NODE" in tip.name:
            recorded_taxa = find_clade_and_move(tree, tip.name)
            if recorded_taxa:
                node_records.append(f'{tip.name}: {"; ".join(recorded_taxa)}')
    with open(node_output_file, 'w') as file:
        for record in node_records:
            file.write(f'{record}\n')
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

def run_command(command, step_name, log_file, critical=False):
    try:
        subprocess.run(command, shell=True, check=True)
        log_status(log_file, f"{step_name}: SUCCESS")
    except subprocess.CalledProcessError:
        log_status(log_file, f"{step_name}: FAILURE")
        print(f"Error: {step_name} failed. Check {log_file} for details.")
        
        if critical:
            exit(1)  # Stop execution if a critical step fails

def process_gene(gene_name_shorter, tree_dir, log_file):
    try:
        tree_list_file = "./{gene_name_shorter}.loop.treelist.txt"
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

def main(threads, proj_name, gene_list, log_file, threshold, use_flag, input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    log_status(log_file, f"Created directory {output_dir}")
    run_command(f'cp {input_dir}/*.tre {output_dir}', f"Copy tree files to {output_dir}", log_file)
    with open(gene_list, 'r') as genes:
        gene_names = [gene.strip() for gene in genes]
    with ThreadPoolExecutor(max_workers=int(threads)) as executor:
        for gene_name_shorter in gene_names:
            executor.submit(process_gene, gene_name_shorter, output_dir, log_file)
    final_output = process_matrices(output_dir, proj_name, threshold, use_flag)
    output_file = f'{output_dir}/{proj_name}.summary_dist.csv'
    final_output.to_csv(output_file, index=False)
    log_status(log_file, f"Processed matrices saved to {output_file}")
    cumulative_output_file = f'{proj_name}.cumulative_dist.csv'
    group_and_sum(output_file, cumulative_output_file)
    log_status(log_file, f"Generated cumulative distance file: {cumulative_output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline for distance calculation and tree processing.")
    parser.add_argument("-t", "--threads", required=True, help="Number of threads to use for parallel processing")
    parser.add_argument("-p", "--proj_name", required=True, help="Project name")
    parser.add_argument("-g", "--gene_list", type=str, default="gene_list.txt", help="Path to the gene list file")
    parser.add_argument("--threshold", type=float, default=1.96, help="Threshold for value adjustment")
    parser.add_argument("--use_flag", action="store_true", help="Use flag to assign the smallest number as 0 and others as 999")
    parser.add_argument("--input_dir", type=str, default="03_phylo_results", help="Input directory containing the tree files")
    parser.add_argument("--output_dir", type=str, default="04_all_trees", help="Output directory for storing results")

    args = parser.parse_args()

    log_file = f"{args.proj_name}_03_distance_calc.log"
    if os.path.exists(log_file):
        os.remove(log_file)
    log_status(log_file, "Pipeline started with the following parameters:")
    log_status(log_file, f"  Threads: {args.threads}")
    log_status(log_file, f"  Project Name: {args.proj_name}")
    log_status(log_file, f"  Gene List: {args.gene_list}")
    log_status(log_file, f"  Threshold: {args.threshold}")
    log_status(log_file, f"  Use Flag: {args.use_flag}")
    log_status(log_file, f"  Input Directory: {args.input_dir}")
    log_status(log_file, f"  Output Directory: {args.output_dir}")

    main(args.threads, args.proj_name, args.gene_list, log_file, args.threshold, args.use_flag, args.input_dir, args.output_dir)

    log_status(log_file, "Pipeline completed successfully.")
    print(f"Pipeline completed. Check {log_file} for the status of each step.")
