
import argparse
import pandas as pd

def process_column(column):
    """Process column data by keeping only the part before the first underscore."""
    return [x.split('_')[0] if '_' in x else x for x in column]

def main():
    # Create an argparse object
    parser = argparse.ArgumentParser(description='Process a file and summarize data.')
    parser.add_argument('input_file', type=str, help='Path to the input file.')
    parser.add_argument('output_file', type=str, help='Path to the output file.')
    args = parser.parse_args()

    # Read the file
    df = pd.read_csv(args.input_file)

    # Process the first column data
    df['processed_col1'] = process_column(df.iloc[:, 0])

    # Summarize the second column data based on the processed first column
    summary = df.groupby('processed_col1')['total_value'].sum().reset_index()

    # Rename the second column of the summary to 'sum_of_second_column' or any other descriptive name
    summary = summary.rename(columns={'total_value': 'total_value'})
    summary = summary.rename(columns={'processed_col1': 'row_name'})

    # Output the summary to a file
    summary.to_csv(args.output_file, index=False)
    print(f'Summary has been written to {args.output_file}')

if __name__ == '__main__':
    main()
