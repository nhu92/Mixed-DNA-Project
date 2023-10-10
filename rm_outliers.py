import pandas as pd
from scipy.stats import f_oneway

# Load the CSV file into a DataFrame
df = pd.read_csv("E:/OneDrive/OneDrive - Texas Tech University/2023Fall/mixed_dna_proj/4471_cb.csv")

# Store the first row and column names for later reference
column_names = df.columns.tolist()
row_names = df.iloc[:, 0].tolist()

# Remove the first column and first row (assuming they contain labels)
df = df.iloc[1:, 1:].apply(pd.to_numeric)

# Perform one-way ANOVA for each column
alpha = 0.001
p_values = []
for col in df.columns:
    group_data = df[col]
    f_statistic, p_value = f_oneway(group_data, df.values.flatten())
    p_values.append(p_value)

# Create a dictionary to store the results
result_dict = {
    "Column": column_names[1:],
    "P-Value": p_values
}
# Create a DataFrame from the results
result_df = pd.DataFrame(result_dict)

# Filter columns with p-values less than alpha to identify statistically different columns
statistically_different_columns = result_df[result_df["P-Value"] < alpha]

print(statistically_different_columns)

# Drop outlier columns from the original DataFrame
df_without_outliers = df.drop(columns=statistically_different_columns["Column"].tolist())

# Save the modified DataFrame to a new CSV file
output_csv_path = "E:/OneDrive/OneDrive - Texas Tech University/2023Fall/mixed_dna_proj/4471_cb_without_outliers.csv"
df_without_outliers.to_csv(output_csv_path, index=False)