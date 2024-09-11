import argparse
import pandas as pd
from scipy.stats import zscore

def process_column(column, taxonomic_level):
    """Process column data by extracting the part according to the taxonomic level (o/f/g/s)."""
    level_map = {'o': 0, 'f': 1, 'g': 2, 's': [2, 3]}
    result = []

    for item in column:
        parts = item.split('_')
        if taxonomic_level == 's':
            result.append('_'.join(parts[2:4]) if len(parts) >= 4 else item)
        else:
            index = level_map[taxonomic_level]
            result.append(parts[index] if len(parts) > index else item)

    return result

def select_taxonomy_by_zscore(summary, zscore_threshold):
    """Select taxonomy names based on a given z-score threshold."""
    return summary.loc[summary['z_score'] > zscore_threshold, 'row_name'].tolist()

def main():
    # Create an argparse object
    parser = argparse.ArgumentParser(description='Process a file and summarize data based on taxonomic level.')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to the output file.')
    parser.add_argument('-tl', '--taxonomic_level', choices=['o', 'f', 'g', 's'], required=True, help='Taxonomic level to process (o: Order, f: Family, g: Genus, s: Species).')
    parser.add_argument('-z', '--zscore_threshold', type=float, required=True, help='Z-score threshold to select taxonomy names.')
    parser.add_argument('-to', '--taxonomy_output_file', type=str, required=True, help='Path to the output file for selected taxonomy names.')
    args = parser.parse_args()

    # Read the file
    df = pd.read_csv(args.input_file)

    # Process the first column data according to the specified taxonomic level
    df['processed_col1'] = process_column(df.iloc[:, 0], args.taxonomic_level)

    # Summarize the second column data based on the processed first column
    summary = df.groupby('processed_col1')['total_value'].sum().reset_index()

    # Calculate the Z-score for the total_value column
    summary['z_score'] = zscore(summary['total_value'])

    # Rename the columns for clarity
    summary = summary.rename(columns={'processed_col1': 'row_name', 'total_value': 'sum_of_total_value'})

    # Sort the result by the sum_of_total_value column
    summary = summary.sort_values(by='sum_of_total_value', ascending=False)

    # Output the summary to a file
    summary.to_csv(args.output_file, index=False)
    print(f'Summary has been written to {args.output_file}')

    # Select taxonomy names based on the z-score threshold
    selected_taxonomies = select_taxonomy_by_zscore(summary, args.zscore_threshold)

    # Write the selected taxonomy names to the output file
    with open(args.taxonomy_output_file, 'w') as f:
        for taxonomy in selected_taxonomies:
            f.write(f"{taxonomy}\n")

    print(f'Selected taxonomy names have been written to {args.taxonomy_output_file}')

if __name__ == '__main__':
    main()
