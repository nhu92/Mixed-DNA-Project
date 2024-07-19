from collections import Counter
import pandas as pd
import numpy as np

def analyze_files(cumulative_csv_path, candidates_tsv_path, output_csv_path):
    # Reading the cumulative CSV file
    df = pd.read_csv(cumulative_csv_path)
    mean = df['total_value'].mean()
    std = df['total_value'].std()

    # Extracting project name and counting individuals
    project_name = cumulative_csv_path.split('/')[-1].split('.cumulative')[0]
    individuals_count = len(project_name.split('x'))

    # Delta calculations
    deltas = np.arange(-1, 2, 0.1)
    results = {}
    for delta in deltas:
        threshold = mean + delta * std
        greater_values = df[df['total_value'] > threshold]['row_name'].tolist()
        results[round(delta, 1)] = greater_values

    # Reading the candidates TSV file
    candidates_df = pd.read_csv(candidates_tsv_path, sep='\t')
    project_ids = [int(id) for id in project_name.split('x')]
    true_related_species = candidates_df[candidates_df['ID'].isin(project_ids)]['related_species'].tolist()
    true_related_species_set = set(true_related_species)

    # Calculating TP, FP, TN, FN, and Rates
    table_results = []
    for delta, row_names in results.items():
        row_names_count = Counter(row_names)
        
        TP = sum(row_names_count[species] for species in true_related_species_set)
        FP = sum(row_names_count[species] for species in row_names_count if species not in true_related_species_set)
        FN = len(true_related_species_set - set(row_names))
        TN = len(set(df['row_name']) - set(row_names) - true_related_species_set)
        PP = len(true_related_species)
        NN = len(set(df['row_name']) - true_related_species_set)

        TPR = TP / (TP + FN) if (TP + FN) > 0 else 0
        FNR = 1 - TPR
        TNR = TN / (TN + FP) if (TN + FP) > 0 else 0
        FPR = 1 - TNR

        # Some other statistics
        PPV = TP / (TP + FP) if (TP + FP) > 0 else 0
        NPV = TN / (FN + TN) if (FN + TN) > 0 else 0
        ACC = (TP + TN) / (PP + NN)

        table_results.append({
            'Individual IDs': 'x'.join([str(id) for id in project_ids]),
            'Individual Number': individuals_count,
            'Delta Value': round(delta, 1),
            'Number of Rows Passed': len(row_names),
            'True Positive Rate': TPR,
            'False Positive Rate': FPR,
            'True Negative Rate': TNR,
            'False Negative Rate': FNR,
            'Positive Predictive Value': PPV,
            'Negative Predictive Value': NPV,
            'Accuracy': ACC
        })

    # Creating a DataFrame from the results and saving it as a CSV file
    results_df = pd.DataFrame(table_results)
    results_df.to_csv(output_csv_path, index=False)
    return "Analysis complete. Results saved to: " + output_csv_path

# Setting up argparse for command line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process cumulative distribution and candidates files.')
    parser.add_argument('cumulative_csv', type=str, help='Path to the cumulative CSV file')
    parser.add_argument('candidates_tsv', type=str, help='Path to the candidates TSV file')
    parser.add_argument('output_csv', type=str, help='Path for the output CSV file')

    args = parser.parse_args()
    result_message = analyze_files(args.cumulative_csv, args.candidates_tsv, args.output_csv)
    print(result_message)
