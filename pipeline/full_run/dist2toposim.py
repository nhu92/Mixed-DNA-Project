import argparse
import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

def distance_to_similarity(distance_df):
    # Copy the DataFrame to avoid changes to the original one
    similarity_df = distance_df.copy()
    # Apply the transformation to convert distance to similarity for numeric columns only
    numeric_cols = similarity_df.select_dtypes(include=[ 'float64', 'int' ]).columns
    similarity_df[numeric_cols] = 1 / (1 + similarity_df[numeric_cols])
    return similarity_df

def clean_up_matrix(df, proj_name, threshold, taxa_file, use_flag=False):

    # Remove rows where row names contain the project name
    df = df[~df[df.columns[0]].str.contains(proj_name)]

    # Keep columns related to the project name
    cols_to_keep = [df.columns[0]] + [col for col in df.columns[1:] if proj_name in col]
    df = df[cols_to_keep]

    # Checking each cell in every column and modify values if needed
    for col in df.columns[1:]:  # Skip the first column as it's often non-numeric (like names, IDs, etc.)
        if use_flag:
            # Assign the smallest number as 0 and the rest as 999
            min_value = df[col].min()
            df[col] = df[col].apply(lambda x: 0 if x == min_value else 999)
        else:
            # Original logic
            col_mean = df[col].mean()
            col_sd = df[col].std()
            df[col] = df[col].apply(lambda x: 999 if x > (col_mean - threshold * col_sd) else x)
        
    if taxa_file:
        # Read the taxa file and store the information in a dictionary
        species_to_taxa = {}
        with open(taxa_file, 'r') as file:
            for line in file:
                parts = line.strip().split(':')
                if len(parts) == 2:
                    species, taxa = parts
                    species_to_taxa[species.strip()] = [tax.strip() for tax in taxa.split(';')]

        # Loop through each pair (column header and taxa list) and update df accordingly
        for header in df.columns[1:]:  # Assuming the first column of df contains the species names
            for species, taxa_list in species_to_taxa.items():
                # Check if the column header matches any species in the taxa list
                if species == header:  # Match the header with the species before ":"
                    for index, row in df.iterrows():  # Loop through rows to check each cell
                        # If the row's species name is not in the taxa list, set its value to 0
                        if row[df.columns[0]] not in taxa_list:  # Compare species names with taxa
                            df.at[index, header] = 999

    
    # Clean up row names by removing digits and trailing underscores
    df.iloc[:, 0] = df.iloc[:, 0].str.replace(r'\d+', '', regex=True).str.rstrip('_')
    return df



def normalize_columns(df):
    # Ensure there are numeric columns to normalize
    if df.select_dtypes(include=[np.number]).empty:
        raise ValueError("No numeric columns to normalize.")

    # Select only the numeric columns for normalization, excluding any non-numeric columns
    numeric_cols = df.select_dtypes(include=[np.number]).columns

    # Initialize the StandardScaler
    scaler = StandardScaler()

    # Normalize the numeric columns
    df[numeric_cols] = scaler.fit_transform(df[numeric_cols])

    return df

def process_matrices(directory, proj_name, threshold, flag):
    all_matrices = []
    
    for filename in os.listdir(directory):
        if filename.endswith('cleaned.csv'):
            # Load the matrix
            matrix_path = os.path.join(directory, filename)
            matrix = pd.read_csv(matrix_path)
            
            # Construct the filename for the corresponding list file
            prefix = filename.split("cleaned.csv")[0]
            list_file = f"{prefix}list.txt"
            list_file_path = os.path.join(directory, list_file)
            # Clean up the matrix
            matrix = clean_up_matrix(matrix, proj_name, threshold, list_file_path, flag) 
            # Apply distance to similarity transformation
            matrix = distance_to_similarity(matrix)

            # Check if there are any numeric columns left before normalization
#             if not matrix.select_dtypes(include=[np.number]).empty:
#                 # Normalize the matrix
#                 matrix = normalize_columns(matrix)

            # Append to the list of all matrices
            all_matrices.append(matrix)

    # Concatenate all matrices into a single DataFrame
    total_matrix = pd.concat(all_matrices, ignore_index=True)

    # Merge all matrices by row names
    total_matrix = pd.concat(all_matrices, ignore_index=True)
    total_matrix.fillna(0, inplace=True)
    # Sum only the numeric columns, excluding the first column with row names
    numeric_cols = [col for col in total_matrix.columns if col != 'Unnamed: 0']
    total_matrix['total_value'] = total_matrix[numeric_cols].sum(axis=1)
    final_output = total_matrix[['Unnamed: 0', 'total_value']].rename(columns={'Unnamed: 0': 'row_name'})
    
    return final_output

def main():
    parser = argparse.ArgumentParser(description='Process batch of matrix files into a single similarity matrix.')
    parser.add_argument('input_dir', type=str, help='Directory containing the matrix CSV files')
    parser.add_argument('output_file', type=str, help='Output file path for the similarity matrix')
    parser.add_argument('removing_pattern', type=str, help='Name of the target sequences to be removed')
    parser.add_argument('--threshold', type=float, default=1.96, help='Threshold for value adjustment in clean_up_matrix function')
    parser.add_argument('--use_flag', action='store_true', help='Use the flag to assign the smallest number as 0 and rest as 999 in clean_up_matrix')
    
    args = parser.parse_args()
    
    final_output = process_matrices(args.input_dir, args.removing_pattern, args.threshold, args.use_flag)
    final_output.to_csv(args.output_file, index=False)
    
    print(f"Processed matrices saved to {args.output_file}")

if __name__ == "__main__":
    main()


