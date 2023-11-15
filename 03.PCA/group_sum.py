import pandas as pd
import argparse

def main(input_file, output_file):
    file_path = input_file
    data = pd.read_csv(file_path)
    
    # Group by the first column and sum the second column
    grouped_sum = data.groupby('row_name')['total_value'].sum().reset_index()
    
    # Save the grouped and summed data to a new CSV file
    grouped_sum.to_csv(output_file, index=False)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a contour plot from a CSV file.")
    parser.add_argument("input_file", help="Input CSV file path")
    parser.add_argument("output_file", help="Output CSV file path")
    args = parser.parse_args()
    main(args.input_file, args.output_file)