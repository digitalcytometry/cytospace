import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import pickle
from scipy.spatial.transform import Rotation
from cytospace.common import read_file

def plot_output(cell_type_path, num_row, num_column, rotation_degrees, rotation_flag, plot_visium, spot_size, plot_marker,
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
                 new_cell_index, index, assigned_nodes, st_path, coords_path, cell_type_path):

    file_delim = "," if coords_path.endswith(".csv") else "\t"
    df_coords = pd.read_csv(coords_path, delimiter=file_delim)
    st_names = list(df_coords.iloc[:,0])
    df_coords.index = st_names

    coords_columns = list(df_coords.columns)
    #if len(coords_columns) > 3:
    #    additional_columns = coords_columns[3:]

    cell_ids_selected_list = cell_ids_selected.tolist()

    assigned_node_names = [st_names[l] for l in index]

    file_delim = "," if cell_type_path.endswith(".csv") else "\t"
    df_labels = pd.read_csv(cell_type_path, delimiter=file_delim)
    df_labels.index = list(df_labels.iloc[:,0])

    numzeros = int(math.log10(len(cell_ids_selected_list)))+1
    unique_cids = ['UCID'+str(i).zfill(numzeros) for i in range(len(cell_ids_selected_list))]

    df_locations = pd.DataFrame.from_dict({'UniqueCID': unique_cids,'OriginalCID': cell_ids_selected_list,
                                'CellType': list(df_labels.loc[cell_ids_selected_list,:][df_labels.columns[1]]),
                                'SpotID': assigned_node_names,
                                'row': list(assigned_locations.iloc[:, 0]),
                                'col': list(assigned_locations.iloc[:, 1])})
    fout = str(output_path)+'/'+str(output_prefix)+'assigned_locations.csv'
    df_locations.to_csv(fout,index=False)
  
    
    metadata = df_locations.copy()
    df = pd.DataFrame(columns=metadata['CellType'].unique(),index=metadata['SpotID'].unique())
    for idx in df.index:
        df.loc[idx] = metadata[metadata['SpotID']==idx]['CellType'].value_counts() 
    df = df.fillna(0)
    df['Total cells'] = df.sum(axis=1)
    df = df.astype(int)
    df.index.name='SpotID'
    fout = str(output_path)+'/'+str(output_prefix)+'cell_type_assignments_by_spot.csv'
    df.to_csv(fout)


    total_cells = np.array(df['Total cells'],dtype=float)
    df = df.iloc[:,:-1].copy()
    df_fracs = df.div(total_cells,axis=0)
    fout = str(output_path)+'/'+str(output_prefix)+'fractional_abundances_by_spot.csv'
    df_fracs.to_csv(fout)
