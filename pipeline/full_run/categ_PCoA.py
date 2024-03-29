import argparse
import csv

def process_row(row, format_type):
    parts = row[0].split('_')
    if format_type == 'order':
        new_first_column = parts[0]
    elif format_type == 'family':
        new_first_column = parts[0][0] + '_' + parts[1]
    elif format_type == 'species':
        new_first_column = parts[2] + '_' + parts[3] if len(parts) >= 4 else 'N/A'
    else:
        raise ValueError("Invalid format type. Choose 'order', 'family', or 'species'.")
    return [new_first_column, row[1], row[2]]

def main():
    parser = argparse.ArgumentParser(description='Process an input file based on the specified format.')
    parser.add_argument('input_file', type=str, help='The input file path.')
    parser.add_argument('output_file', type=str, help='The output file path.')
    parser.add_argument('format_type', type=str, choices=['order', 'family', 'species'], help='The format to apply to the first column.')
    
    args = parser.parse_args()

    with open(args.input_file, 'r') as infile, open(args.output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        
        # Handle header if present
        header = next(reader)
        writer.writerow(header)  # Write the header to the output file
        
        for row in reader:
            processed_row = process_row(row, args.format_type)
            writer.writerow(processed_row)

if __name__ == "__main__":
    main()
