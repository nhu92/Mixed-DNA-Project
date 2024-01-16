import pandas as pd
import argparse

def process_files(cumulative_dist_path, candidates_path, output_path):
    # Read the CSV and TSV files
    cumulative_dist_df = pd.read_csv(cumulative_dist_path)
    candidates_df = pd.read_csv(candidates_path, sep='\t')

    # Extract individual IDs from the cumulative distribution filename
    extracted_ids = cumulative_dist_path.split('/')[-1].split('.')[0].split('x')
    extracted_ids = [int(id_) for id_ in extracted_ids]  # Convert to integers

    # Calculate the mean - 0.2 * standard deviation
    mean_value = cumulative_dist_df['total_value'].mean()
    std_dev = cumulative_dist_df['total_value'].std()
    threshold = mean_value - 0.2 * std_dev

    # Find the rows where total_value is greater than this threshold
    outliers = cumulative_dist_df[cumulative_dist_df['total_value'] > threshold]
    outlier_names = outliers['row_name'].tolist()

    # Match these names with the individual IDs in the candidates file
    results = []
    for id_ in extracted_ids:
        related_species = candidates_df[candidates_df['ID'] == id_]['related_species'].iloc[0]
        status = "successful" if related_species in outlier_names else "failed"
        results.append({"ID": id_, "Status": status})

    # Convert results to a DataFrame and save to output file
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, index=False)
    print(f"Processed data saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description='Process files and output the result.')
    parser.add_argument('cumulative_dist_path', type=str, help='Path to the cumulative distribution file')
    parser.add_argument('candidates_path', type=str, help='Path to the candidates file')
    parser.add_argument('output_path', type=str, help='Path for the output file')

    args = parser.parse_args()

    process_files(args.cumulative_dist_path, args.candidates_path, args.output_path)

if __name__ == "__main__":
    main()