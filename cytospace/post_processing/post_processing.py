import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
from scipy.spatial.transform import Rotation
from cytospace.common import read_file


def save_results(output_path, output_prefix, cell_ids_selected, cell_ids_new, assigned_locations,
                 new_cell_index, index, assigned_nodes, st_path, coords_path,
                 cell_type_path, assigned_locations_path):

    file_delim = "," if coords_path.endswith(".csv") else "\t"
    df_coords = pd.read_csv(coords_path, delimiter=file_delim)
    st_names = list(df_coords.iloc[:,0])
    df_coords.index = st_names

    coords_columns = list(df_coords.columns)
    if len(coords_columns) > 3:
        additional_columns = coords_columns[3:]

    cell_ids_selected_list = cell_ids_selected.tolist()
    cell_ids_new_list = cell_ids_new.tolist()
    assigned_node_names = [st_names[l] for l in index]

    file_delim = "," if cell_type_path.endswith(".csv") else "\t"
    df_labels = pd.read_csv(cell_type_path, delimiter=file_delim)
    df_labels.index = list(df_labels.iloc[:,0])

    numzeros = int(math.log10(len(cell_ids_selected_list)))+1
    unique_cids = ['UCID'+str(i).zfill(numzeros) for i in range(len(cell_ids_selected_list))]

    df_locations = pd.DataFrame.from_dict({'UniqueCID': unique_cids,'OriginalCID': cell_ids_selected_list,
                                'PlaceHolderCID': cell_ids_new_list,
                                'CellType': list(df_labels.loc[cell_ids_selected_list,:][df_labels.columns[1]]),
                                'SpotID': assigned_node_names,
                                'row': list(assigned_locations.iloc[:, 0]),
                                'col': list(assigned_locations.iloc[:, 1])})
    df_locations.to_csv(assigned_locations_path, index=False)
    
    metadata = df_locations.copy()
    df = pd.DataFrame(columns=metadata['CellType'].unique(),index=metadata['SpotID'].unique())
    for idx in df.index:
        df.loc[idx] = metadata[metadata['SpotID']==idx]['CellType'].value_counts() 
    df = df.fillna(0)
    df['Total cells'] = df.sum(axis=1)
    df = df.astype(int)
    df.index.name='SpotID'
    fout = output_path / f'{output_prefix}cell_type_assignments_by_spot.csv'
    df.to_csv(fout)


    total_cells = np.array(df['Total cells'],dtype=float)
    df = df.iloc[:,:-1].copy()
    df_fracs = df.div(total_cells,axis=0)
    fout = output_path / f'{output_prefix}fractional_abundances_by_spot.csv'
    df_fracs.to_csv(fout)

