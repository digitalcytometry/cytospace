import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import pickle
from scipy.spatial.transform import Rotation
from cytospace.common import read_file


def plot_output(cell_type_path, num_row, num_column, rotation_degrees, rotation_flag, spot_size,
                output_path, output_prefix, assigned_nodes, new_cell_index, coordinates):
    number_of_cells_per_page = num_row * num_column
    cell_type_labels = read_file(cell_type_path)
    cell_type_labels_unique = ['All cells'] + list(cell_type_labels.columns)
    clusters = len(cell_type_labels_unique)
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
                node_assignment[i] = np.sum(assignment_cell_type == i)
            x = coordinates[:, 0]
            y = coordinates[:, 1]
            if max(node_assignment) == 0:
               ps = plt.scatter(x, y, s=155, c=node_assignment, vmin=0, vmax=5, marker='h')
            else:
               ps = plt.scatter(x, y, s=155, c=node_assignment, marker='h')    
            len_x = (np.max(x) - np.min(x))
            len_y = (np.max(y) - np.min(y))
            x_limit = (127 - len_x)/2
            y_limit = (77 - len_y)/2
            plt.xlim(np.min(x) - x_limit, np.max(x) + x_limit)
            plt.ylim(np.min(y) - y_limit, np.max(y) + y_limit)         
            plt.title(cell_type_labels_unique[j + number_of_cells_per_page*s], fontsize=80)
            plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Arial']})
            plt.tight_layout()

            cb = plt.colorbar(ps, shrink=0.85)
            cb.ax.tick_params(labelsize=60)
            ax.axis("off")

            index = s + 1    
            plt.savefig(output_path / f"{output_prefix}plot_cell_type_locations_{index}.pdf")

        clusters = clusters - number_of_cells_per_page
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
