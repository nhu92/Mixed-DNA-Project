import argparse
import pandas as pd
import numpy as np
from scipy import stats


def analyze_files(cumulative_csv_path, candidates_tsv_path, output_csv_path):
    # Reading the cumulative CSV file
    df = pd.read_csv(cumulative_csv_path)
    mean = df['total_value'].mean()
    std = df['total_value'].std()

    # Extracting project name and counting individuals
    project_name = cumulative_csv_path.split('/')[-1].split('.cumulative')[0]
    individuals_count = len(project_name.split('x'))

    
    original_data = df['total_value']
    
    # Alpha values to test
    alpha_values = [0.0001, 0.001, 0.005, 0.01, 0.05, 0.1]
    
    # Dictionary to store the final results
    results_item = {}
    
    # Analysis for each alpha value
    for alpha in alpha_values:
        data_sorted = df.sort_values(by='total_value', ascending=False)
        row_names = []  # List to store removed rows as outliers
    
        # Iteratively remove the highest 'total_value' and perform one-tailed t-test
        for i in range(len(data_sorted)):
            if not data_sorted.empty:
                removed_row = data_sorted.iloc[0]
                data_sorted = data_sorted.iloc[1:]
                row_names.append(removed_row['row_name'])

    
                # Perform one-tailed t-test
                if len(data_sorted) > 1:
                    t_stat, p_value = stats.ttest_ind(original_data, data_sorted['total_value'], equal_var=False)
                    one_tailed_p_value = p_value / 2  # One-tailed p-value
                else:
                    one_tailed_p_value = np.nan  # Not enough data
    
                # Check for significant difference
                if one_tailed_p_value < alpha:
                    break

        # Add the results to the final dictionary
        results_item[alpha] = row_names
        
    # Reading the candidates TSV file
    candidates_df = pd.read_csv(candidates_tsv_path, sep='\t')
    project_ids = [int(id) for id in project_name.split('x')]
    true_related_species = candidates_df[candidates_df['ID'].isin(project_ids)]['related_species'].tolist()

    # Calculating TP, FP, TN, FN, and Rates
    table_results = []
    for alpha, row_names in results_item.items():
        TP = len(set(row_names).intersection(set(true_related_species)))
        print(row_names, true_related_species)
        FP = len(set(row_names) - set(true_related_species))
        TN = len(set(df['row_name']) - set(row_names) - set(true_related_species))
        FN = len(set(true_related_species) - set(row_names))
        PP = len(true_related_species)
        NN = len(set(df['row_name']) - set(true_related_species))

        TPR = TP / len(true_related_species) if true_related_species else 0
        FNR = 1 - TPR
        TNR = TN / NN
        FPR = 1 - TNR

        # Some other statistics
        PPV = TP / (TP + FP) if (TP + FP) > 0 else 0
        NPV = TN / (FN + TN) if (FN + TN) > 0 else 0
        ACC = (TP + TN) / (PP + NN)


        table_results.append({
            'Individual IDs': 'x'.join([str(id) for id in project_ids]),
            'Individual Number': individuals_count,
            'Delta Value': round(alpha, 4),
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
