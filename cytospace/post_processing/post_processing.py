import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import pickle
from scipy.spatial.transform import Rotation
from cytospace.common import read_file


def plot_output(cell_type_path, num_row, num_column, rotation_degrees, rotation_flag,
                output_path, output_prefix, assigned_nodes, new_cell_index, coordinates):
    cell_type_labels = read_file(cell_type_path)
    cell_type_labels_unique = cell_type_labels.columns
    clusters = cell_type_labels_unique.shape[0]
    num_rows = int(round(np.sqrt(clusters)))
    plt.figure(figsize=(num_rows * 2.5 * 3.5, num_rows * 2.5 * 3))
    noSpots = coordinates.shape[0]

    rotation_radians = np.radians(rotation_degrees)
    rotation_axis = np.array([0, 0, 1])
    rotation_vector = rotation_radians * rotation_axis
    rotation = Rotation.from_rotvec(rotation_vector)
    rotated_coordinates = np.zeros((coordinates.shape[0], 2))
    for k in range(coordinates.shape[0]):
        rotated_vec = rotation.apply([coordinates.iloc[k, 0], coordinates.iloc[k, 1], 0])
        rotated_coordinates[k, :] = rotated_vec[0:2]

    coordinates = rotated_coordinates if rotation_flag else coordinates.values
    for counter, j in enumerate(range(clusters)):
        ax = plt.subplot(num_row, num_column, counter + 1)

        # Calculate number of assigned cells to each spot
        assignment_cell_type = assigned_nodes[int(new_cell_index[j]):int(new_cell_index[j + 1] - 1)]
        node_assignment = np.zeros(noSpots)
        for i in range(noSpots):
            node_assignment[i] = np.sum(assignment_cell_type == i)

        ps = plt.scatter(coordinates[:, 0], coordinates[:, 1],
                         s=(240 / clusters) * int(clusters / 4), c=node_assignment, marker='H')
        plt.xlim(np.min(coordinates[:, 0]), np.max(coordinates[:, 0]))
        plt.ylim(np.min(coordinates[:, 1]), np.max(coordinates[:, 1]))

        plt.title(cell_type_labels_unique[j], fontsize=25)
        plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Arial']})
        plt.tight_layout()

        cb = plt.colorbar(ps, shrink=0.71)
        cb.ax.tick_params(labelsize=25)
        ax.axis("off")

        plt.savefig(output_path / f"{output_prefix}plot_cell_type_locations.pdf")

    plt.close()


def save_results(output_path, output_prefix, cell_ids_selected, assigned_locations,
                 new_cell_index, index, assigned_nodes):
    # Store python objects of the final result for plotting or further processing
    for result, name in [(assigned_nodes, "assigned_nodes"), (index, "index"),
                         (new_cell_index, "new_cell_index"),
                         (assigned_locations, "assigned_locations")]:
        with open(output_path / f"{output_prefix}{name}.pckl", "wb") as f:
            pickle.dump(result, f)

    # Save assigned locations of the up/down sampled cells
    cell_ids_selected_list = cell_ids_selected.tolist()
    index_df = pd.DataFrame(index)
    with open(output_path / f"{output_prefix}assigned_nodes.csv", 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(('Cell IDs', 'Node'))
        for row in range(cell_ids_selected.shape[0]):
            writer.writerow([cell_ids_selected_list[row], index_df.iloc[row,0] + 1])

    with open(output_path / f"{output_prefix}assigned_locations.csv", 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(('Cell IDs', 'X', 'Y'))
        for row in range(cell_ids_selected.shape[0]):
            writer.writerow([cell_ids_selected_list[row],
                             assigned_locations.iloc[row, 0],
                             assigned_locations.iloc[row, 1]])

    with open(output_path / f"{output_prefix}new_cell_type.csv", 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow('n')
        for row in range(new_cell_index.shape[0]):
            writer.writerow([new_cell_index[row]])
