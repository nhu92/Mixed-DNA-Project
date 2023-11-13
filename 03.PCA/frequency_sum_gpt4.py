import argparse
import pandas as pd
import matplotlib.pyplot as plt

def read_csv(file_path, header='infer'):
    return pd.read_csv(file_path, header=header)

def plot_graphs(data, color_mapping, taxa_mapping, output_bar, output_pie):
    # Counting the frequency of items in the second column
    item_counts = data.iloc[:, 1].value_counts()

    # Map colors and taxa to items
    colors = item_counts.index.map(lambda x: color_mapping.get(x, 'gray'))
    taxa = item_counts.index.map(lambda x: taxa_mapping.get(x, 'Unknown'))

    # Create a legend mapping from taxa to color
    unique_taxa = set(taxa)
    legend_elements = [plt.Line2D([0], [0], color=color_mapping[t], lw=4, label=t) for t in unique_taxa if t in color_mapping]

    # Plotting Bar Graph
    plt.figure(figsize=(15, 10))  # Increased size
    item_counts.plot(kind='bar', color=colors)
    plt.title('Item Frequency Bar Graph')
    plt.xlabel('Items')
    plt.ylabel('Frequency')
    plt.legend(handles=legend_elements, title='Taxa')
    plt.savefig(output_bar, format='svg')

    # Plotting Pie Chart
    plt.figure(figsize=(15, 10))  # Increased size
    item_counts.plot(kind='pie', colors=colors, autopct='%1.1f%%')
    plt.title('Item Frequency Pie Chart')
    plt.legend(handles=legend_elements, title='Taxa', bbox_to_anchor=(1, 1))  # Adjusted position for legend
    plt.savefig(output_pie, format='svg')

def main():
    parser = argparse.ArgumentParser(description='Plot item frequency graphs from CSV data.')
    parser.add_argument('data_file', help='CSV file without header containing the data')
    parser.add_argument('color_file', help='CSV file with item color mappings')
    parser.add_argument('output_bar', help='Output path for the bar graph SVG file')
    parser.add_argument('output_pie', help='Output path for the pie chart SVG file')
    args = parser.parse_args()

    # Read data and color mapping files
    data = read_csv(args.data_file, header=None)
    color_data = read_csv(args.color_file)

    # Check if the color_data DataFrame has the expected columns
    expected_columns = ['item', 'color', 'taxa']
    if not all(col in color_data.columns for col in expected_columns):
        raise ValueError(f"Color mapping file must contain columns: {expected_columns}")

    # Create color and taxa mapping dictionaries
    color_mapping = dict(zip(color_data['item'], color_data['color']))
    taxa_mapping = dict(zip(color_data['item'], color_data['taxa']))

    # Plot and save graphs
    plot_graphs(data, color_mapping, taxa_mapping, args.output_bar, args.output_pie)

if __name__ == '__main__':
    main()
