import argparse
import pandas as pd
import re

def read_files(csv_file, tsv_file):
    # Read CSV and TSV files
    csv_data = pd.read_csv(csv_file)
    tsv_data = pd.read_csv(tsv_file, sep='\t')
    return csv_data, tsv_data

def extract_ids_from_filename(filename):
    # Extract IDs from filename using regular expression
    return re.findall(r'\d+', filename)

def calculate_threshold(csv_data):
    # Calculate mean and standard deviation, then determine the threshold
    mean_value = csv_data['total_value'].mean()
    std_dev_value = csv_data['total_value'].std()
    return mean_value - 0.2 * std_dev_value

def identify_outliers(csv_data, threshold):
    # Identify outliers in the CSV data
    return csv_data[csv_data['total_value'] > threshold]

def match_ids_and_species(csv_outliers, tsv_data):
    # Match outliers with related species in the TSV data
    outlier_species = csv_outliers['row_name'].tolist()
    return tsv_data[tsv_data['related_species'].isin(outlier_species)]

def determine_identification_status(file_ids, matched_ids):
    # Determine identification status for each ID
    success_ids = set(matched_ids['ID'])
    failed_ids = set(file_ids) - success_ids
    success_results = pd.DataFrame({'ID': list(success_ids), 
                                    'identification_status': 'Successfully Identified'})
    fail_results = pd.DataFrame({'ID': list(failed_ids), 
                                 'identification_status': 'Fail to be Identified'})
    return pd.concat([success_results, fail_results])

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process CSV and TSV files to identify species.')
    parser.add_argument('csv_file', type=str, help='Input CSV file')
    parser.add_argument('tsv_file', type=str, help='Input TSV file')
    parser.add_argument('output_file', type=str, help='Output file for results')

    # Parse arguments
    args = parser.parse_args()

    # Read input files
    csv_data, tsv_data = read_files(args.csv_file, args.tsv_file)

    # Extract IDs from CSV filename
    file_ids = extract_ids_from_filename(args.csv_file)
    file_ids = [int(id) for id in file_ids]

    # Calculate threshold and identify outliers
    threshold = calculate_threshold(csv_data)
    csv_outliers = identify_outliers(csv_data, threshold)

    # Match IDs and species, then determine identification status
    matched_ids = match_ids_and_species(csv_outliers, tsv_data)
    results = determine_identification_status(file_ids, matched_ids)

    # Write results to the output file
    results.to_csv(args.output_file, index=False)

if __name__ == '__main__':
    main()
