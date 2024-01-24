import argparse
import pandas as pd

def calculate_delta_ppv(cumulative_dist_path, delta_ppv_path, output_path):
    # Load the data
    cumulative_dist_df = pd.read_csv(cumulative_dist_path)
    delta_ppv_df = pd.read_csv(delta_ppv_path)

    # Calculate mean and standard deviation
    mean = cumulative_dist_df['total_value'].mean()
    std_dev = cumulative_dist_df['total_value'].std()

    # Calculate delta
    cumulative_dist_df['delta'] = (cumulative_dist_df['total_value'] - mean) / std_dev

    # Function to find the closest PPV
    def find_closest_ppv(delta):
        delta_rounded = round(delta, 2)
        closest_match = delta_ppv_df.iloc[(delta_ppv_df['Delta Value'] - delta_rounded).abs().argsort()[:1]]
        return closest_match['Positive Predictive Value'].values[0]

    # Apply the function to find PPV
    cumulative_dist_df['ppv'] = cumulative_dist_df['delta'].apply(find_closest_ppv)

    # Save to CSV
    cumulative_dist_df.to_csv(output_path, index=False)
    print(f"Output file saved as {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Process cumulative distribution and delta PPV data.")
    parser.add_argument("cumulative_dist_file", help="Path to the cumulative distribution CSV file")
    parser.add_argument("delta_ppv_file", help="Path to the delta PPV CSV file")
    parser.add_argument("output_file", help="Path for the output CSV file")

    args = parser.parse_args()

    calculate_delta_ppv(args.cumulative_dist_file, args.delta_ppv_file, args.output_file)

if __name__ == "__main__":
    main()
