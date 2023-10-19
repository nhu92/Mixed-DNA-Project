import argparse
import pandas as pd

# Function to calculate mean and standard deviation within a node
def calculate_node_stats(data, node):
    subset = data.filter(like=node, axis=1)
    mean = subset.mean(axis=1)
    std = subset.std(axis=1)
    return mean, std

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Calculate mean and standard deviation of genetic distances within each node.")
    parser.add_argument("input_file", help="Input file containing genetic distances")
    parser.add_argument("output_file", help="Output file for storing calculated statistics")
    args = parser.parse_args()

    # Read input data from the CSV file
    data = pd.read_csv(args.input_file, index_col=0)

    # Rename columns by removing ".{number}" from the column names
    data.columns = data.columns.str.replace(r'\.\d+', '', regex=True)

    # Get unique node names
    unique_nodes = set(data.columns[1:])

    # Create a new DataFrame for the output
    output_data = pd.DataFrame(index=data.index)

    for node in unique_nodes:
        mean, std = calculate_node_stats(data, node)
        output_data[f"{node}_mean"] = mean
        output_data[f"{node}_sd"] = std

    # Save the output data to a CSV file
    output_data.to_csv(args.output_file)

if __name__ == "__main__":
    main()
