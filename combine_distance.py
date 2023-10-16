import pandas as pd
import numpy as np
import glob
import re
import argparse

def remove_rows(df, remove_string):
    # Remove rows containing the input string in the first column and empty first column
    df = df[~df[df.columns[0]].str.contains(remove_string, case=False, na=False, regex=True)]
    df = df[df[df.columns[0]].notna()]
    return df.reset_index(drop=True)

def process_tsv_files(input_pattern, output_file, remove_string):
    # Step 1: Read all TSV files matching the input pattern
    tsv_files = glob.glob(input_pattern)
    
    if not tsv_files:
        print("No matching TSV files found.")
        return

    dfs = []
    for file in tsv_files:
        df = pd.read_csv(file, delimiter=',')
        df = remove_rows(df, remove_string)
        df = df.sort_values(by=df.columns[0])  # Sort by the first column
        dfs.append(df)

    # Step 2: Row bind all input TSV files
    merged_df = pd.concat(dfs, axis=1, ignore_index=False)
    # Step 3: Remove columns named "Unnamed: 0"

    first_column = merged_df.iloc[:, 0]
    merged_df = merged_df.loc[:, ~merged_df.columns.str.contains('^Unnamed: 0$', case=False)]

    merged_df.index = first_column
    merged_df.to_csv(f"{output_file}.tsv", sep='\t', index=True)
	
    # Step 4: Group columns by the pattern "*_NODE{number}_*"
    node_pattern = re.compile(r'_(NODE_\d+)_')
    new_column_names = {'species': first_column}
    
    for col in merged_df.columns:
        col_name = str(col)
        match = node_pattern.search(col_name)
        if match:
            node_number = match.group(1)
            new_column_names[col_name] = f"NODE{node_number}"

    merged_df = merged_df.rename(columns=new_column_names)
    merged_df.to_csv(f"{output_file}.tsv", sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sort, group, calculate mean and standard deviation, and output results to a TSV file.")
    parser.add_argument("--input-pattern", required=True, help="Input pattern (wildcard) to match TSV files.")
    parser.add_argument("--output-file", required=True, help="Output TSV file for results.")
    parser.add_argument("--remove-string", default="", help="String to remove from the first column.")
    args = parser.parse_args()

    input_pattern = args.input_pattern
    output_file = args.output_file
    remove_string = args.remove_string

    process_tsv_files(input_pattern, output_file, remove_string)
