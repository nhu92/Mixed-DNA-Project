import re
import argparse

def remove_pattern_from_newick(newick_str):
    # Remove all instances of '_R_'
    newick_str = re.sub(r'_R_', '', newick_str)
    
    # Remove the last '_number' pattern
    newick_str = re.sub(r'_\d+', '', newick_str)

    return newick_str

def main():
    parser = argparse.ArgumentParser(description='Remove specified patterns from a Newick phylogeny file')
    parser.add_argument('input_file', help='Input Newick phylogeny file')
    parser.add_argument('output_file', help='Output file for cleaned phylogeny')

    args = parser.parse_args()

    # Read the Newick phylogeny from the input file
    with open(args.input_file, 'r') as file:
        newick_data = file.read()

    # Remove the specified patterns
    newick_data_cleaned = remove_pattern_from_newick(newick_data)

    # Write the cleaned data to the output file
    with open(args.output_file, 'w') as file:
        file.write(newick_data_cleaned)

    print(f'Cleaned data saved to {args.output_file}')

if __name__ == '__main__':
    main()
