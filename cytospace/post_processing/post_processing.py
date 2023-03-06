import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
from scipy.spatial.transform import Rotation
from cytospace.common import read_file

def save_results(output_path, output_prefix, cell_ids_selected, all_cells_save, assigned_locations,
                 cell_type_data, assigned_locations_path, sampling_method, single_cell):
    """
    Parameters :
        output_path, output_prefix  (str) : /path/to/output/dir and outprefix_ parts, respectively.
        assigned_locations_path     (str) : /path/to/save/assigned_locations.csv

        cell_ids_selected  (1D np.array[str], List[str], or other list-like type) 
        assigned_locations (pd.DataFrame) :
            cell_ids_selected is a list of single cell IDs (from scRNA_data.columns).
            assigned_locations is a subset of coordinates_data. index=SpotID; columns=row/col, X/Y, etc.
            the orders of these two have to match, in that the nth cell_id in cell_ids_selected
                is assigned to the spot corresponding to the nth row in assigned_locations.

        all_cells_save     (pd.DataFrame)
        save_sc_expression         (bool) :
            all_cells_save is the unnormalized expression matrix of sampled single cells,
                including the original and placeholder (if applicable). index=CellID; columns=CellType.
            this matrix is saved as a CSV file if save_sc_expression is True.
        
        cell_type_data     (pd.DataFrame) : index=CellID; columns=celltype. as read in from read_data().
        sampling_method             (str) : as specified in args.
    
    Returns None.
    """

    cell_ids_selected_list = list(cell_ids_selected)

    numzeros = int(math.log10(len(cell_ids_selected_list)))+1
    unique_cids = ['UCID'+str(i).zfill(numzeros) for i in range(len(cell_ids_selected_list))]
    
    # if a cell is from the original sc dataset : original_cell_ids = single-cell cell ID, cell_types = corresponding cell type
    # if a cell is generated as a placeholder   : original_cell_ids = 'CELL_NA', cell_types = retrieved from placeholder cell ID
    original_cell_ids = [cell_id
                            if (cell_id in cell_type_data.index)
                            else 'CELL_NA'
                            for cell_id in cell_ids_selected_list]
    cell_types        = [cell_type_data.loc[cell_id, :][cell_type_data.columns[0]]
                            if cell_id in cell_type_data.index
                            else f"TYPE_{cell_id.split('_')[1]}"
                            for cell_id in cell_ids_selected_list]
    
    cell_ids_selected_list = [cid[5:] for cid in cell_ids_selected_list]    # remove 'CELL_'
    assigned_node_names    = [sid[5:] for sid in assigned_locations.index]  # remove 'SPOT_'
    original_cell_ids      = [cid[5:] for cid in original_cell_ids]         # remove 'CELL_'
    cell_types             = [tid[5:] for tid in cell_types]                # remove 'TYPE_'

    if sampling_method == "place_holders":
        
        df_locations = pd.DataFrame.from_dict({'UniqueCID': unique_cids,'OriginalCID': original_cell_ids,
                                    'PlaceHolderCID': cell_ids_selected_list,
                                    'CellType': cell_types,
                                    'SpotID': assigned_node_names,
                                    assigned_locations.columns.values[0]: list(assigned_locations.iloc[:, 0]),
                                    assigned_locations.columns.values[1]: list(assigned_locations.iloc[:, 1])})
    else:
        
        df_locations = pd.DataFrame.from_dict({'UniqueCID': unique_cids,'OriginalCID': original_cell_ids,
                                    'CellType': cell_types,
                                    'SpotID': assigned_node_names,
                                    assigned_locations.columns.values[0]: list(assigned_locations.iloc[:, 0]),
                                    assigned_locations.columns.values[1]: list(assigned_locations.iloc[:, 1])})
    
    df_locations.to_csv(assigned_locations_path, index=False)
    if sampling_method == "place_holders":
        fout_scrna = output_path / f'{output_prefix}new_scRNA.csv'
        all_cells_save = all_cells_save.loc[:, ~all_cells_save.columns.isin(cell_type_data.index)]
        all_cells_save.index   = [gene_id[5:] for gene_id in all_cells_save.index.astype(str)]   # remove 'GENE_'
        all_cells_save.columns = [cell_id[5:] for cell_id in all_cells_save.columns.astype(str)] # remove 'CELL_'
        all_cells_save.to_csv(fout_scrna)

    if not single_cell:
        # count each pair of SpotID x CellType, then expand to a 2D count matrix
        df = df_locations.loc[:, ['SpotID', 'CellType']].value_counts().unstack(fill_value=0)\
                .reindex(index=df_locations.SpotID.unique(), columns=df_locations.CellType.unique())
        df['Total cells'] = df.sum(axis=1)
        df = df.astype(int)
        df.index.name='SpotID'
        fout = output_path / f'{output_prefix}cell_type_assignments_by_spot.csv'
        df.to_csv(fout)

        total_cells = np.array(df['Total cells'],dtype=float)
        df = df.iloc[:, :-1]
        df_fracs = df.div(total_cells,axis=0)
        fout = output_path / f'{output_prefix}fractional_abundances_by_spot.csv'
        df_fracs.to_csv(fout)
