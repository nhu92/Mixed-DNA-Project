import subprocess
import os
import argparse
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor

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

                # Run matrix_ult.py
                matrix_cmd = f'python matrix_ult.py -t "{filename}" -n "{tree_dir}/{gene_name_shorter}.{i}.list.txt" -o "{tree_dir}/{gene_name_shorter}.{i}.matrix"'
                run_command(matrix_cmd, f"Run matrix_ult.py for {gene_name_shorter} tree {i}", log_file)

                # Copy matrix to cleaned CSV
                copy_cmd = f'cp "{tree_dir}/{gene_name_shorter}.{i}.matrix" "{tree_dir}/{gene_name_shorter}.{i}.cleaned.csv"'
                run_command(copy_cmd, f"Copy matrix to cleaned CSV for {gene_name_shorter} tree {i}", log_file)

        os.remove(tree_list_file)
        log_status(log_file, "Removed temporary file loop.treelist.txt")
    except Exception as e:
        log_status(log_file, f"Failed processing {gene_name_shorter}: {e}")
        print(f"Failed processing {gene_name_shorter}: {e}")

def main(threads, proj_name, gene_list, log_file):
    tree_dir = "04_all_trees"
    os.makedirs(tree_dir, exist_ok=True)
    log_status(log_file, "Created directory ./all_trees")

    run_command(f'cp 03_phylo_results/*.tre {tree_dir}', "Copy tree files to ./04_all_trees", log_file)

    # Process each gene in parallel
    with open(gene_list, 'r') as genes:
        gene_names = [gene.strip() for gene in genes]

    with ThreadPoolExecutor(max_workers=int(threads)) as executor:
        for gene_name_shorter in gene_names:
            executor.submit(process_gene, gene_name_shorter, tree_dir, log_file)

    # Run dist2toposim.py
    dist2toposim_cmd = f'python dist2toposim.py {tree_dir} ./{proj_name}.summary_dist.csv {proj_name} --threshold 1'
    run_command(dist2toposim_cmd, "Run dist2toposim.py", log_file)

    # Run group_sum.py
    group_sum_cmd = f'python group_sum.py ./{proj_name}.summary_dist.csv ./{proj_name}.cumulative_dist.csv'
    run_command(group_sum_cmd, "Run group_sum.py", log_file)

    # Generate Z-axis file
    grep_cmd = f'grep -v {proj_name} {proj_name}.cumulative_dist.csv | cut -d \',\' -f2 > {proj_name}.Zaxis.csv'
    run_command(grep_cmd, "Generate Z-axis file", log_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline for distance calculation and tree processing.")
    parser.add_argument("-t", "--threads", required=True, help="Number of threads to use for parallel processing")
    parser.add_argument("-p", "--proj_name", required=True, help="Project name")
    parser.add_argument("-g", "--gene_list", required=True, help="Path to the gene list file")

    args = parser.parse_args()

    log_file = f"{args.proj_name}_03_distance_calc.log"

    log_status(log_file, "Pipeline started with the following parameters:")
    log_status(log_file, f"  Threads: {args.threads}")
    log_status(log_file, f"  Project Name: {args.proj_name}")
    log_status(log_file, f"  Gene List: {args.gene_list}")

    main(args.threads, args.proj_name, args.gene_list, log_file)

    log_status(log_file, "Pipeline completed successfully.")
    print(f"Pipeline completed. Check {log_file} for the status of each step.")
