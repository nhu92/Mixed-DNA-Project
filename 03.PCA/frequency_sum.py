import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Create bar and pie charts from a CSV file")
    parser.add_argument("input_file", help="Input CSV file")
    parser.add_argument("output_svg", help="Output SVG file")
    args = parser.parse_args()

    # Read the CSV file into a pandas DataFrame without a header
    df = pd.read_csv(args.input_file, header=None)

    # Get the second column
    column = df.iloc[:, 1]

    # Count the frequency of each item in the second column
    item_counts = column.value_counts()

    # Create a colorful bar chart
    colors = cm.viridis(np.linspace(0, 1, len(item_counts)))
    item_counts.plot(kind="bar", color=colors)
    plt.title("Bar Chart of Second Column")
    plt.xlabel("Item")
    plt.ylabel("Frequency")
    plt.savefig("bar_chart.svg", format="svg")
    plt.clf()  # Clear the current figure

    # Create a colorful pie chart with a legend
    labels = item_counts.index
    sizes = item_counts.values
    explode = [0.1] * len(labels)  # Adjust the value (0.1) for desired effect
    plt.figure()
    plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct="%1.1f%%", startangle=140)
    plt.axis("equal")
    plt.title("Pie Chart of Second Column")

    # Create a legend for the pie chart on the right side
    legend_labels = [f"{label} ({size})" for label, size in zip(labels, sizes)]
    legend = plt.legend(legend_labels, title="Items", bbox_to_anchor=(1, 0.5), loc="center left")
    plt.savefig(args.output_svg, format="svg", bbox_extra_artists=(legend,))

if __name__ == "__main__":
    main()
