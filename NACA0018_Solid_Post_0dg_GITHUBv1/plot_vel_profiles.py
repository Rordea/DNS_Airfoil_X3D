# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 22:59:07 2024

@author: ricar
"""
import os
import matplotlib.pyplot as plt
import xcompact3d_toolbox as x3d
from matplotlib import rcParams


rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.size'] = 18


filename_baseline = './results_data/vel_profiles.dat' 
compare_filename = '../NACA0018_Prs_Post_0dg/results_data/vel_profiles.dat'


def read_and_plot_profiles(baseline_path, compare_path):
    def read_data(file_path):
        """Reads and transposes data from the given file path."""
        headers = []
        data = []

        with open(file_path, mode='r') as file:
            # Read the header line
            headers = file.readline().strip().split('\t\t')
            
            # Read the data lines
            for line in file:
                # Split the line into individual values
                values = list(map(float, line.strip().split()))
                data.append(values)

        # Transpose the data to get columns
        data_transposed = list(map(list, zip(*data)))
        return headers, data_transposed

    # Read data from both datasets
    baseline_headers, baseline_data = read_data(baseline_path)
    compare_headers, compare_data = read_data(compare_path)

    # Extract y-coordinates (first column)
    y_coords = baseline_data[0]

    # Set up the figure for subplots
    num_profiles = len(baseline_headers) - 1
    rows = 1
    cols = 6

    # Ensure there are enough subplots; if not, adjust layout
    total_subplots = min(num_profiles, rows * cols)
    fig, axes = plt.subplots(rows, cols, figsize=(15, 10), sharey=True)

    # Flatten the axes array for easy iteration
    axes = axes.flatten()

    # Plot each velocity profile on a separate subplot
    for idx, profile_name in enumerate(baseline_headers[1:], start=1):
        if idx <= total_subplots:
            axes[idx - 1].plot(
                baseline_data[idx], y_coords, label=f"Solid X/C={profile_name}",
                color='darkred', marker='v', linestyle='--', markevery=25
            )
            axes[idx - 1].plot(
                compare_data[idx], y_coords, label=f"Porous X/C={profile_name}",
                color='black', marker='x', linestyle='--', markevery=25
            )
            #xes[idx - 1].set_xlabel('Velocity')
            #axes[idx - 1].set_ylabel('Y coordinates')
            #axes[idx - 1].set_title(profile_name)
            axes[idx - 1].grid(True)
            #axes[idx - 1].legend(loc='best')

    # Hide any unused subplots
    for ax in axes[num_profiles:]:
        ax.set_visible(False)

    # Adjust layout for the subplot grid
    #plt.tight_layout()

    # Create a new figure for all profiles together
    plt.figure(figsize=(10, 8))
    for idx, profile_name in enumerate(baseline_headers[1:], start=1):
        plt.plot(
            baseline_data[idx], y_coords, label=f"Solid X/C={profile_name}",
            color='darkred', marker='v', linestyle='--', markevery=25
        )
        plt.plot(
            compare_data[idx], y_coords, label=f"Porous X/C={profile_name}",
            color='blue', marker='o', linestyle='-', markevery=25
        )

    # Labeling and displaying the combined plot


# Call the function with the file path
read_and_plot_profiles(filename_baseline,compare_filename)


plt.xlabel('Velocity u')
plt.ylabel('Y coordinates')
plt.title('All Velocity Profiles Combined')
plt.grid(True)
#plt.legend(loc='best')
plt.minorticks_on()
plt.grid(which='both', color='gray', linestyle='-', linewidth=0.5, alpha=0.2)
#plt.legend(framealpha=0.0)
#legend = plt.legend( loc='best', ncol=2,fontsize=25,frameon=False)
plt.tight_layout()
plt.show()
#filename_to_compare = os.path.join(os.path.abspath(os.path.join(os.getcwd(), os.pardir)), f'NACA0018_Prs_Post_{AoA}dg/results_data/forces.dat' )