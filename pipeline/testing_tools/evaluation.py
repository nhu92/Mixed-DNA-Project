# Function to perform the analysis with two input files and save the output as a CSV file
def analyze_files(cumulative_csv_path, candidates_tsv_path, output_csv_path):
    # Reading the cumulative CSV file
    df = pd.read_csv(cumulative_csv_path)
    mean = df['total_value'].mean()
    std = df['total_value'].std()

    # Extracting project name and counting individuals
    project_name = cumulative_csv_path.split('/')[-1].split('.cumulative')[0]
    individuals_count = len(project_name.split('x'))

    # Delta calculations
    deltas = np.arange(-1, 1.1, 0.1)
    results = {}
    for delta in deltas:
        threshold = mean + delta * std
        greater_values = df[df['total_value'] > threshold]['row_name'].tolist()
        results[round(delta, 1)] = greater_values

    # Reading the candidates TSV file
    candidates_df = pd.read_csv(candidates_tsv_path, sep='\t')
    project_ids = [int(id) for id in project_name.split('x')]
    true_related_species = candidates_df[candidates_df['ID'].isin(project_ids)]['related_species'].tolist()

    # Calculating TP, FP, TN, FN, and Rates
    table_results = []
    for delta, row_names in results.items():
        TP = len(set(row_names).intersection(set(true_related_species)))
        FP = len(set(row_names) - set(true_related_species))
        TN = len(set(df['row_name']) - set(row_names) - set(true_related_species))
        FN = len(set(true_related_species) - set(row_names))

        FPR = FP / (FP + TN) if (FP + TN) > 0 else 0
        TNR = TN / (TN + FP) if (TN + FP) > 0 else 0
        FNR = FN / (FN + TP) if (FN + TP) > 0 else 0
        TPR = TP / len(true_related_species) if true_related_species else 0

        table_results.append({
            'Individual Number': individuals_count,
            'Delta Value': round(delta, 1),
            'Number of Rows Passed': len(row_names),
            'True Positive Rate': TPR,
            'False Positive Rate': FPR,
            'True Negative Rate': TNR,
            'False Negative Rate': FNR
        })

    # Creating a DataFrame from the results and saving it as a CSV file
    results_df = pd.DataFrame(table_results)
    results_df.to_csv(output_csv_path, index=False)
    return "Analysis complete. Results saved to: " + output_csv_path

# Example usage of the function
# analyze_files('/path/to/cumulative_dist.csv', '/path/to/candidates.tsv', '/path/to/output.csv')
# Note: Replace '/path/to/...' with actual file paths when using this function.
