import argparse
import pandas as pd
from sklearn import metrics
from sklearn.preprocessing import LabelEncoder
import re

def calculate_statistics(result_file, true_file, output_file):
    # Load data
    result = pd.read_csv(result_file)
    result = result[~result['Sample'].str.startswith('NODE')]
    result['Sample'] = result['Sample'].apply(lambda x: re.sub(r'_[0-9]{4}_mean', '', x))
    true = pd.read_csv(true_file)


    # Merge the result and true dataframes on 'Sample'
    merged = pd.merge(result, true, on='Sample', how='left')

    # Use label encoder to ensure that the cluster labels match
    le = LabelEncoder()
    merged['Cluster_x'] = le.fit_transform(merged['Cluster_x'])
    merged['Cluster_y'] = le.fit_transform(merged['Cluster_y'])

    # Calculate statistics
    ari = metrics.adjusted_rand_score(merged['Cluster_y'], merged['Cluster_x'])
    ami = metrics.adjusted_mutual_info_score(merged['Cluster_y'], merged['Cluster_x'])
    v_measure = metrics.v_measure_score(merged['Cluster_y'], merged['Cluster_x'])

    # Save to output file
    with open(output_file, 'w') as f:
        f.write(f"{result_file}, {ari}, {ami}, {v_measure}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Evaluate KMeans clustering results.')
    parser.add_argument('--result_file', type=str, help='The file containing the clustering results.')
    parser.add_argument('--true_file', type=str, help='The file containing the true clustering.')
    parser.add_argument('--output_file', type=str, help='The output file to save the statistics.')
    
    args = parser.parse_args()
    
    calculate_statistics(args.result_file, args.true_file, args.output_file)
