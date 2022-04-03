
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
from scipy.spatial.transform import Rotation
import matplotlib.pyplot as plt

from cytospace.common import read_file, add_plotting_arguments


def argument_parser():
    parser = argparse.ArgumentParser(description='Plot CytoSPACE results')

    # Required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument("-alp", "--assigned-locations-path", type=str, required=True,
                          help="Path to the assigned locations by CytoSPACE", default=None)
    required.add_argument("-cp", "--coordinates-path", type=str, required=True,
                         help="Path to transcriptomics data (coordinates)", default=None)
    parser.add_argument("-o", "--output-filename", type=str, required=True,
                        help="Output file name, i.e. 'cytospace_results.pdf'",  default=None)

    # Plotting
    add_plotting_arguments(parser)

    arguments = parser.parse_args()

    return arguments.__dict__


def plot_results(assigned_locations_path, coordinates_path, output_filename, num_row=4, num_column=4, rotation_flag=True, plot_visium=True,
                 rotation_degrees=270, spot_size=175, plot_marker='h'):
    # Read data
    coordinates = read_file(coordinates_path)
    coordinates_id = coordinates.index.values
    assigned_locations = read_file(assigned_locations_path)
    assigned_nodes = assigned_locations.loc[:,'SpotID']

    number_of_cells_per_page = num_row * num_column
    cell_type_labels = assigned_locations.loc[:,'CellType']
    cell_type_labels_unique = ['All cells'] + list(sorted(np.unique(cell_type_labels), key=str.lower))
    clusters = len(cell_type_labels_unique)
    new_cell_index = np.zeros((clusters))
    sum_temp = 0
    for l in range(1,clusters):
        sum_temp = sum_temp + np.sum(cell_type_labels == cell_type_labels_unique[l])
        new_cell_index[l] = sum_temp
        
    noSpots = coordinates.shape[0]
    rotation_radians = np.radians(rotation_degrees)
    rotation_axis = np.array([0, 0, 1])
    rotation_vector = rotation_radians * rotation_axis
    rotation = Rotation.from_rotvec(rotation_vector)
    rotated_coordinates = np.zeros((coordinates.shape[0], 2))
    for k in range(coordinates.shape[0]):
        rotated_vec = rotation.apply([coordinates.iloc[k, 0], coordinates.iloc[k, 1], 0])
        rotated_coordinates[k, :] = rotated_vec[0:2]

    # Set 'global' plot options
    plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Arial']})

    coordinates = rotated_coordinates if rotation_flag else coordinates.values
    new_cell_index = list(new_cell_index)
    new_cell_index.insert(0, 0)
    for s in range(int(clusters/number_of_cells_per_page) + 1):
        plt.figure(figsize=(52, 52))
        for counter, j in enumerate(range(min(number_of_cells_per_page, clusters))):
            ax = plt.subplot(num_row, num_column, counter + 1)

            # Calculate number of assigned cells to each spot
            if s == 0 and counter == 0:
                assignment_cell_type = assigned_nodes
            else:
                assignment_cell_type = assigned_nodes[int(new_cell_index[j + number_of_cells_per_page*s]):int(new_cell_index[j + number_of_cells_per_page*s + 1] - 1)]
            node_assignment = np.zeros(noSpots)
            for i in range(noSpots):
                node_assignment[i] = np.sum(assignment_cell_type == coordinates_id[i])
            x = coordinates[:, 0]
            y = coordinates[:, 1]
            if max(node_assignment) == 0:
               ps = plt.scatter(x, y, s=spot_size, c=node_assignment, vmin=0, vmax=5, marker=plot_marker)
            else:
               ps = plt.scatter(x, y, s=spot_size, c=node_assignment, marker=plot_marker)

            if plot_visium:
                len_x = (np.max(x) - np.min(x))
                len_y = (np.max(y) - np.min(y))
                x_limit = (127 - len_x)/2
                y_limit = (77 - len_y)/2
            else:
                len_x = (np.max(x) - np.min(x))
                len_y = (np.max(y) - np.min(y))
                x_limit = len_x/10
                y_limit = len_y/10
            plt.xlim(np.min(x) - x_limit, np.max(x) + x_limit)
            plt.ylim(np.min(y) - y_limit, np.max(y) + y_limit)
            plt.title(cell_type_labels_unique[j + number_of_cells_per_page*s], fontsize=80)
            plt.tight_layout()

            cb = plt.colorbar(ps, shrink=0.85)
            cb.ax.tick_params(labelsize=60)
            ax.axis("off")

            index = s + 1
            output_filename = Path(output_filename)
            plt.savefig(output_filename.parent / (output_filename.stem + f'_{index}' + output_filename.suffix))

        clusters = clusters - number_of_cells_per_page
        plt.close()


def run_cytospace_plot():
    arguments = argument_parser()
    plot_results(**arguments)


if __name__ == "__main__":
    arguments = argument_parser()
    plot_results(**arguments)