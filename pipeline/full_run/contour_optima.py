import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter, maximum_filter, label
from scipy.spatial import KDTree
import argparse

def load_data(file_path):
    return pd.read_csv(file_path)

def find_local_maxima(x, y, z):
    tree = KDTree(list(zip(x, y)))
    radius = (max(x) - min(x)) / 20
    maxima_points = []
    for i in range(len(z)):
        indices = tree.query_ball_point([x[i], y[i]], r=radius)
        if all(z[i] >= z[j] for j in indices) and z[i] >= 0:
            maxima_points.append((x[i], y[i]))
    return maxima_points

def create_contour_plot(data, output_file):
    # Extracting X, Y, and Z values
    x = data.iloc[:, 1].values
    y = data.iloc[:, 2].values
    z = data.iloc[:, 3].values
    names = data.iloc[:, 0].values

    # Grid for interpolation
    grid_x, grid_y = np.mgrid[min(x):max(x):1000j, min(y):max(y):1000j]
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='linear')
    zi_smoothed = gaussian_filter(grid_z, sigma=1)

    # Identifying local maxima
    local_maxima_points = find_local_maxima(x, y, z)

    # Plotting
    plt.figure(figsize=(12, 10))
    contour = plt.contourf(grid_x, grid_y, zi_smoothed, levels=15, cmap='terrain')
    plt.colorbar(contour)
    plt.scatter(x, y, marker='o', c='white', s=10)
    for i, txt in enumerate(names):
        plt.annotate(txt, (x[i], y[i]), fontsize=8, ha='right')
    for x_peak, y_peak in local_maxima_points:
        plt.scatter(x_peak, y_peak, color='red', s=50)
        index = np.where((x == x_peak) & (y == y_peak))[0][0]
        plt.annotate(names[index], (x_peak, y_peak), color='red', fontsize=12, ha='right')
    plt.title('Contour Plot with Positive Local Maxima and Annotations')
    plt.xlabel('X axis')
    plt.ylabel('Y axis')
    plt.savefig(output_file, format='svg')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Generate a contour plot from a CSV file.')
    parser.add_argument('input_file', type=str, help='Input CSV file path')
    parser.add_argument('output_file', type=str, help='Output SVG file path')
    args = parser.parse_args()

    data = load_data(args.input_file)
    create_contour_plot(data, args.output_file)

if __name__ == '__main__':
    main()
