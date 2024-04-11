import pandas as pd
import numpy as np
from pandas.core.common import flatten

import random
import time
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


import concurrent.futures
import os
import math

from cytospace.common import read_file, read_visium, normalize_data, downsample, check_paths, argument_parser, estimate_cell_type_fractions
from cytospace.post_processing import save_results, plot_results
from cytospace.linear_assignment_solvers import (calculate_cost, match_solution, import_solver,
                                                 call_solver)

def read_data(scRNA_path, cell_type_path, cell_type_fraction_estimation_path, n_cells_per_spot_path, 
                st_cell_type_path, output_path, output_prefix, spaceranger_path=None, st_path=None, coordinates_path=None):
    if spaceranger_path is not None:
        st_data, coordinates_data = read_visium(spaceranger_path,output_path)
    elif (st_path is None) and (coordinates_path is None):
        raise ValueError("For ST data, you must provide either a tar.gz file or paths for expression and coordinates.")
    else:
        # Read data
        st_data = read_file(st_path)
        coordinates_data = read_file(coordinates_path)

    st_data = st_data[~st_data.index.duplicated(keep=False)]
    if (st_cell_type_path is None) and (cell_type_fraction_estimation_path is None):
        print('Estimating cell type fractions')
        if (spaceranger_path is not None):
            st_data_outpath = os.path.join(output_path, f"{output_prefix}ST_expression.txt")
            st_data.to_csv(st_data_outpath, sep='\t')
            cell_type_fraction_estimation_path = estimate_cell_type_fractions(scRNA_path, cell_type_path, st_data_outpath, output_path, output_prefix)
        else:
            cell_type_fraction_estimation_path = estimate_cell_type_fractions(scRNA_path, cell_type_path, st_path, output_path, output_prefix)

    st_data.columns = ['SPOT_'+str(col) for col in st_data.columns]
    st_data.index = ['GENE_'+str(idx) for idx in st_data.index]

    coordinates_data.index = ['SPOT_'+str(idx) for idx in coordinates_data.index]

    scRNA_data = read_file(scRNA_path)
    scRNA_data.columns = ['CELL_'+str(col) for col in scRNA_data.columns]
    scRNA_data.index = ['GENE_'+str(idx) for idx in scRNA_data.index]

    scRNA_data = scRNA_data[~scRNA_data.index.duplicated(keep=False)]

    cell_type_data = read_file(cell_type_path)
    cell_type_data.index = ['CELL_'+str(idx) for idx in cell_type_data.index]
    cell_type_data.iloc[:,0] = ['TYPE_'+str(cell) for cell in list(cell_type_data.iloc[:,0])]

    #st_data = st_data.loc[(st_data!=0).any(1), (st_data!=0).any(0)]
    #scRNA_data = scRNA_data.loc[(scRNA_data!=0).any(1), (scRNA_data!=0).any(0)]

    if st_cell_type_path is not None:
        st_cell_type_data = read_file(st_cell_type_path)
        st_cell_type_data.index = ['SPOT_'+str(idx) for idx in st_cell_type_data.index]
        st_cell_type_data.iloc[:,0] = ['TYPE_'+str(cell) for cell in list(st_cell_type_data.iloc[:,0])]
    else:
        st_cell_type_data = None
    if cell_type_fraction_estimation_path is not None:
        cell_type_fraction_data = read_file(cell_type_fraction_estimation_path)
        cell_type_fraction_data.columns = ['TYPE_'+str(col) for col in cell_type_fraction_data.columns]
    else:
        cell_type_fraction_data = None
    if n_cells_per_spot_path is not None:
        n_cells_per_spot_data = read_file(n_cells_per_spot_path)
        n_cells_per_spot_data.index = ['SPOT_'+str(idx) for idx in n_cells_per_spot_data.index]
    else:
        n_cells_per_spot_data = None

    # Order data to match
    try:
        st_data = st_data[coordinates_data.index]
        scRNA_data = scRNA_data[cell_type_data.index]
        if st_cell_type_data is not None:
            st_cell_type_data = st_cell_type_data.loc[coordinates_data.index, :]
        if n_cells_per_spot_data is not None:
            n_cells_per_spot_data = n_cells_per_spot_data.transpose(copy=False)
            n_cells_per_spot_data = n_cells_per_spot_data[coordinates_data.index]
            n_cells_per_spot_data = n_cells_per_spot_data.transpose(copy=False)
    except Exception as e:
        raise IndexError(f"The ST data: {st_path} and coordinates data: {coordinates_path} have to "
                         "have the same spot IDs for columns and rows, respectively, "
                         f"and scRNA data: {scRNA_path} and cell type data: {cell_type_path} have"
                         " to have the same cell IDs for columns and rows, respectively.")

    # Validate input
    if (st_data.columns != coordinates_data.index).any():
        raise IndexError(f"The ST data: {st_path} and coordinates data: {coordinates_path} have to "
                         "have the same spot IDs for columns and rows, respectively.")

    if (scRNA_data.columns != cell_type_data.index).any():
        raise IndexError(f"The scRNA data: {scRNA_path} and cell type data: {cell_type_path} have"
                         " to have the same cell IDs for columns and rows, respectively.")
    
    if (st_cell_type_data is not None) and (st_cell_type_data.index != coordinates_data.index).any():
        raise IndexError(f"The ST cell type data: {st_cell_type_path} and coordinates data: {coordinates_path} have to "
                         "have the same spot IDs for rows.")

    if (st_cell_type_data is None) and (cell_type_fraction_data is None):
        raise ValueError("At least one of st_cell_type_path and cell_type_fraction_estimation_path should be specified."
                         "For --single-cell, st_cell_type_path is recommended; if not --single-cell, cell_type_fraction_estimation_path is required.")
    
    if (st_cell_type_data is not None) and (cell_type_fraction_data is not None):
        print("Warning: st_cell_type_path and cell_type_fraction_estimation_path are both specified.")
        print("If --single-cell, cell_type_fraction_estimation_path will be ignored in this case.")

    return scRNA_data, cell_type_data, st_data, coordinates_data, cell_type_fraction_data, n_cells_per_spot_data, st_cell_type_data


def estimate_cell_number_RNA_reads(st_data, mean_cell_numbers):
    # Read data
    expressions = st_data.values.astype(float)

    # Data normalization
    expressions_tpm_log = normalize_data(expressions)

    # Set up fitting problem
    RNA_reads = np.sum(expressions_tpm_log, axis=0, dtype=float)
    mean_RNA_reads = np.mean(RNA_reads)
    min_RNA_reads = np.min(RNA_reads)

    min_cell_numbers = 1 if min_RNA_reads > 0 else 0

    fit_parameters = np.polyfit(np.array([min_RNA_reads, mean_RNA_reads]),
                                np.array([min_cell_numbers, mean_cell_numbers]), 1)
    polynomial = np.poly1d(fit_parameters)
    cell_number_to_node_assignment = polynomial(RNA_reads).astype(int)

    return cell_number_to_node_assignment


def get_cell_type_fraction(number_of_cells, cell_type_fraction_data):
    # Uncomment commented lines for closer numbers
    cell_type_numbers = cell_type_fraction_data.transpose()
    #fractions_to_numbers = np.round(cell_type_numbers.values*number_of_cells)
    fractions_to_numbers = cell_type_numbers.values*number_of_cells
    fractions_to_numbers = fractions_to_numbers.astype(int)
    cell_type_numbers.iloc[:,0] = fractions_to_numbers
    max_ct = cell_type_numbers.idxmax().values[0]
    #cell_type_numbers.loc[max_ct,cell_type_numbers.columns[0]] += number_of_cells - sum(cell_type_numbers.iloc[:,0])
    cell_type_numbers.loc[cell_type_numbers.index[0],cell_type_numbers.columns[0]] += number_of_cells - sum(cell_type_numbers.iloc[:,0])
    return cell_type_numbers

def partition_indices(indices, split_by_category_list=None, split_by_interval_int=None, shuffle=True):
    """
    Splits the provided indices into list of smaller index sets based on other parameters.
    indices is originally a single 1D numpy array, which is then split and returned as a list of smaller 1D numpy arrays.
    e.g., indices = np.arange(0, 5000) can be split into [np.arange(0, 2000), np.arange(2000, 5000)].

    Parameters :
        indices                (1D np.array(int)) : indices to be split.

        split_by_category_list (1D np.array(int)) : number of indices for each category that cannot be mixed together.
            split_by_category_list should sum to len(indices)
            e.g., if split_by_category_list is [3000, 5000], with indices == 0:8000,
                    then the first 3000 will be partitioned separately from the latter 5000.
                    i.e., [0:2000, 2000:3000, 3000:7000, 7000:8000] is possible, but [0:2000, 2000:7000, 7000:8000] is not.
        
        split_by_interval_int               (int) : max length of each partition.
            if split_by_category_list is None, then indices are split into partitions of this size.
            e.g., if split_by_interval_int is 1000, with indices == 0:2500, then [0:1000, 1000:2000, 2000:2500] is returned.
            if split_by_category_list is specified, then any category that exceeds this size will be further partitioned.
            e.g., if split_by_category_list is [500, 1000, 300] and split_by_interval_int is 400, with indices == 0:1800,
                    then [0:400, 400:500, 500:900, 900:1300, 1300:1500, 1500:1800] is returned.
                    - split_by_category_list sets breakpoints at 500 and 1500
                    - split_by_interval_int sets breakpoints at every 400 inside each of the three groups

        shuffle                            (bool) : whether indices should be shuffled before being split.
    
    Returns :
        List[1D np.array(int)] : Partitioned indices.
    """
    num_indices = len(indices)

    if shuffle:
        np.random.shuffle(indices)

    # initialize breakpoints, as start and end
    breakpoints = [0, num_indices]

    # add breakpoints set by split_by_category_list
    if split_by_category_list is not None:
        if np.sum(split_by_category_list) != num_indices:
            print('Warning: sum of counts in each category does not match the full length')
        breakpoints.extend(np.cumsum(split_by_category_list))
        breakpoints = sorted(np.unique(breakpoints))
    
    # add breakpoints set by split_by_interval_int
    if split_by_interval_int is not None:
        breakpoints_new = breakpoints.copy()
        for idx in range(len(breakpoints)-1):
            if breakpoints[idx+1] - breakpoints[idx] <= split_by_interval_int:
                continue
            breakpoints_new.extend(np.arange(breakpoints[idx], breakpoints[idx+1], split_by_interval_int))
        breakpoints = breakpoints_new
    
    # take unique values, sort in ascending order, and remove start and end indices
    breakpoints = sorted(np.unique(breakpoints))[1:-1]

    # partition the indices based on breakpoints
    split_indices_list = np.array_split(indices, breakpoints)

    return split_indices_list


def sample_single_cells(scRNA_data, cell_type_data, cell_type_numbers_int, sampling_method, seed):
    """
    Samples cells from scRNA_data based on the cell type distribution specified in cell_type_numbers_int.
    The sampled count for each cell type will match the number specified in cell_type_numbers_int.

    Parameters :
        scRNA_data             (2D pd.DataFrame) : gene x cell scRNA expression data to be sampled from.
        cell_type_data        (nx1 pd.DataFrame) : cell types corresponding to scRNA_data. (index=CellID)
        cell_type_numbers_int (nx1 pd.DataFrame) : ST cell count for each cell type. (index=CellType)
        sampling_method (str), seed (int) : as specified in args.
    
    Returns :
        all_cells_save         (2D pd.DataFrame) : gene x cell scRNA expression data for sampled single cells.
            This will be a subset (likely with duplicate columns) of the provided scRNA_data (sampling_method == "duplicates"),
            otherwise a superset of the provided scRNA_data with newly generated placeholder cells (sampling_method == "place_holders").
            It is guaranteed that the cells will be in order of cell types as specified in cell_type_numbers_int.index.
            index=gene, columns=CellID.
    """    
    np.random.seed(seed)
    random.seed(seed)

    # Down/up sample of scRNA-seq data according to estimated cell type fractions
    # follow the order of cell types in cell_type_numbers_int
    unique_cell_type_labels = cell_type_numbers_int.index.values

    # initialize variables
    all_cells_save_list = [] # List of 2D np.ndarray of single cell expression
    cell_names_list = [] # List of 1D np.array of single cell IDs
    sampled_index_total = [] # List of 1D np.array of single cell indices

    for cell_type in unique_cell_type_labels:
        cell_type_index = np.nonzero(cell_type_data.values[:, 0] == cell_type)[0].tolist()
        cell_type_count_available = len(cell_type_index)
        if cell_type_count_available == 0:
            raise ValueError(f"Cell type {cell_type} in the ST dataset is not available in the scRNA-seq dataset.")
        cell_type_count_desired = cell_type_numbers_int.loc[cell_type][0]

        if sampling_method == "place_holders":
            if cell_type_count_desired > cell_type_count_available:
                num_genes = scRNA_data.shape[0]
                num_placeholder_cells = cell_type_count_desired - cell_type_count_available
                scRNA_original_np = scRNA_data.iloc[:, cell_type_index].to_numpy()
                scRNA_placeholder_np = np.zeros((num_genes, num_placeholder_cells))

                sampled_index = np.random.choice(cell_type_index, num_placeholder_cells)

                for i1 in range(num_placeholder_cells):
                    scRNA_placeholder_np[:, i1] = [np.random.choice(scRNA_original_np[j1, :]) for j1 in range(num_genes)]

                # # alternate implementation (vectorization on one axis)
                # for gene_idx in range(num_genes):
                #     rand_idxs = np.random.randint(0, cell_type_count_available, size=num_placeholder_cells)
                #     scRNA_placeholder_np[gene_idx, :] = scRNA_original_np[gene_idx, rand_idxs.tolist()]
                
                cell_names_list.append(np.array(scRNA_data.columns.values[cell_type_index]))
                cell_names_list.append(np.array([cell_type.replace('TYPE_', 'CELL_') + '_new_' + str(i+1) for i in range(num_placeholder_cells)]))
                all_cells_save_list.append(scRNA_original_np)
                all_cells_save_list.append(scRNA_placeholder_np)

            else:
                cell_type_selected_index = random.sample(cell_type_index, cell_type_count_desired)

                cell_names_list.append(scRNA_data.columns.values[cell_type_selected_index])
                all_cells_save_list.append(scRNA_data.iloc[:, cell_type_selected_index].to_numpy())

        elif sampling_method == "duplicates":
            if cell_type_count_desired > cell_type_count_available:
                cell_type_selected_index = np.concatenate([
                    cell_type_index, np.random.choice(cell_type_index, cell_type_count_desired - cell_type_count_available)
                ], axis=0) # ensure at least one copy of each, then sample the rest

            else:
                cell_type_selected_index = random.sample(cell_type_index, cell_type_count_desired)
        
            sampled_index_total.append(cell_type_selected_index)
        
        else:
            raise ValueError("Invalid sampling_method provided")
    
    if sampling_method == "place_holders":
        all_cells_save = pd.DataFrame(
                            np.concatenate(all_cells_save_list, axis=1),
                            index=scRNA_data.index,
                            columns=np.concatenate(cell_names_list, axis=0)
        )
    else:
        sampled_index_total = np.concatenate(sampled_index_total, axis=0).astype(int)
        all_cells_save = scRNA_data.iloc[:, sampled_index_total]

    return all_cells_save


def solve_linear_assignment_problem(scRNA_norm_data, st_norm_data, cell_number_to_node_assignment,
                                    solver_method, solver, seed, distance_metric, process_idx=None):
    """
    Parameters :
        scRNA_norm_data (2D np.ndarray(float)) : normalized gene x cell scRNA expression data to be used as reference
        st_norm_data    (2D np.ndarray(float)) : normalized gene x cell ST expression data to be used as target
        cell_number_to_node_assignment (1D np.ndarray(int)) :
            estimated cell count at each ST spot, where the value at index n denotes the cell count at the nth spot of st_norm_data
        solver_method, solver, seed, distance_metric : as specified in the input args
        process_idx           (int) : returned as is; used to keep track of asynchronous processes in apply_linear_assignment
    Returns :
        mapped_st_index (List[int]) : indices of ST spots (column index of st_norm_data) where each single cell is mapped.
                                        list has length scRNA_norm_data.shape[1]; order of single cells follows scRNA_norm_data.
        process_idx           (int) : returned as is.
    """
    distance_repeat, location_repeat =\
        calculate_cost(scRNA_norm_data, st_norm_data, cell_number_to_node_assignment,
                       solver_method, distance_metric)

    if solver_method == 'lapjv' or solver_method == 'lapjv_compat':
        print('Solving linear assignment problem ...')
        np.random.seed(seed)
        cost_scaled = distance_repeat + 1e-16 * np.random.rand(distance_repeat.shape[0],
                                                               distance_repeat.shape[1])
        t0 = time.perf_counter()
        assignment = call_solver(solver, solver_method, cost_scaled)
        print(f"Time to solve linear assignment problem: {round(time.perf_counter() - t0, 2)} seconds")
        assigned_nodes = location_repeat[assignment]
        mapped_st_index = np.transpose(assigned_nodes).tolist()

    elif solver_method == 'lap_CSPR':
        print('Solving linear assignment problem ...')
        np.random.seed(seed)
        cost_scaled = 10**6 * distance_repeat + 10 * np.random.rand(distance_repeat.shape[0],
                                                                    distance_repeat.shape[1]) + 1
        cost_scaled = np.transpose(cost_scaled)
        cost_scaled_int = cost_scaled.astype(int)
        cost_scaled_int_list = cost_scaled_int.tolist()
        t0 = time.perf_counter()
        assignment = match_solution(cost_scaled_int_list)

        print(f"Time to solve linear assignment problem: {round(time.perf_counter() - t0, 2)} seconds")
        assigned_nodes = location_repeat[assignment[:, 0].astype(int)]
        mapped_st_index = assigned_nodes.tolist()
    else:
        raise ValueError("Invalid solver_method provided")

    return mapped_st_index, process_idx


def apply_linear_assignment(scRNA_data, st_data, coordinates_data, cell_number_to_node_assignment,
                                solver_method, solver, seed, distance_metric, number_of_processors,
                                index_sc_list, index_st_list=None, subsampled_cell_number_to_node_assignment_list=None):
    """
    Parallelizes process by queueing a subprocess for each subset.
    Each subset is specified by index_st_list or subsampled_cell_number_to_node_assignment_list (in reference to st_data)
        along with index_sc_list (in reference to sc_data).
    The output of each subprocess is aggregated and returned as a single object.

    Parameters :
        scRNA_data, st_data, coordinates_data (pd.DataFrame) :
            formatted as read in from read_data().
            scRNA_data will have been sampled upstream.
        cell_number_to_node_assignment (1D np.ndarray(int)) :
            estimated cell count at each ST spot, where the value at index n denotes the cell count at the nth spot of st_data
        solver_method, solver, seed, distance_metric, number_of_processors :
            as specified in the input args.

        index_sc_list (List[1D np.array(int)]) : 
            partition of scRNA_data cell indices (each of length number_of_selected_(sub_)spots, with the exception of the last),
            that denotes the subsets that scRNA_data will be partitioned into for parallelization.

        index_st_list (List[1D np.array(int)]) :
            used if --single-cell.
            partition of st_data spot indices (each of length number_of_selected_spots, with the exception of the last partition),
            that denotes the subsets that st_data will be partitioned into for parallelization.
        subsampled_cell_number_to_node_assignment_list (List[1D np.ndarray(int)]) :
            used if --sampling-sub-spots.
            list of cell_number_to_node_assignment (each of length [spot count of ST data] and
                summing to number_of_selected_sub_spots with the exception of the last),
            where each is used for a subprocess mapping (sampled) scRNA_data to st_data.

    Returns :
        assigned_locations (pd.DataFrame) :
            list of ST coordinates (from coordinates_data) where each single cell in cell_ids_selected is mapped to.
            index = Spot ID; columns = coordinates_data.columns
        cell_ids_selected (1D np.array) :
            list of single cell IDs (from scRNA_data.index; with duplicates if necessary) that are mapped.
            the nth single cell in cell_ids_selected is mapped to the spot on the nth row of assigned_locations.
    """
    if (index_st_list is not None) and (subsampled_cell_number_to_node_assignment_list is not None):
        raise ValueError("index_st_list and subsampled_cell_number_to_node_assignment_list cannot both be specified")
    
    # normalize data; output is an np.ndarray
    scRNA_norm_np = normalize_data(scRNA_data.to_numpy())
    st_norm_np = normalize_data(st_data.to_numpy())

    # regenerate pandas dataframe from the normalized data
    scRNA_norm_data = pd.DataFrame(scRNA_norm_np, index=scRNA_data.index, columns=scRNA_data.columns)
    st_norm_data = pd.DataFrame(st_norm_np, index=st_data.index, columns=st_data.columns)

    if (index_st_list is None) and (subsampled_cell_number_to_node_assignment_list is None):
        mapped_st_index, _ =\
            solve_linear_assignment_problem(
                scRNA_norm_data.iloc[:, index_sc_list[0]].to_numpy(), st_norm_data.to_numpy(), cell_number_to_node_assignment,
                    solver_method, solver, seed, distance_metric)
        assigned_locations = coordinates_data.iloc[mapped_st_index]
        cell_ids_selected = scRNA_norm_data.columns.values[index_sc_list[0]]

        return assigned_locations, cell_ids_selected

    results = []
    cell_ids_selected_list = [] # List of 1D np.array (single cell ID) from each process
    assigned_locations_list = [] # List of pd.DataFrame (ST spot coordinates) from each process

    # compute the number of processes
    if index_st_list is not None:
        # called for --single_cell
        num_iters = len(index_st_list)
    elif subsampled_cell_number_to_node_assignment_list is not None:
        # called for --sampling-sub-spots
        num_iters = len(subsampled_cell_number_to_node_assignment_list)
    else:
        raise ValueError("Invalid point")
    print(f"Number of required processors: {num_iters}")

    with concurrent.futures.ProcessPoolExecutor(max_workers=min(num_iters, number_of_processors)) as executor:
        for idx in range(num_iters):
            if index_st_list is not None:
                # called for --single-cell
                st_norm_data_selected = st_norm_data.iloc[:, index_st_list[idx]]
                cell_number_to_node_assignment_selected = cell_number_to_node_assignment[index_st_list[idx]]
            elif subsampled_cell_number_to_node_assignment_list is not None:
                # called for --sampling-sub-spots
                st_norm_data_selected = st_norm_data
                cell_number_to_node_assignment_selected = subsampled_cell_number_to_node_assignment_list[idx]
            else:
                raise ValueError("Invalid point")

            scRNA_norm_data_selected = scRNA_norm_data.iloc[:, index_sc_list[idx]]
            
            # launch process
            results.append(executor.submit(
                solve_linear_assignment_problem,
                    scRNA_norm_data_selected.to_numpy(), st_norm_data_selected.to_numpy(), cell_number_to_node_assignment_selected,
                    solver_method, solver, seed, distance_metric, process_idx=idx
                )
            )
        
        for f in concurrent.futures.as_completed(results):
            # aggregate results
            mapped_st_index, process_idx = f.result()

            # because the nth single cell in cell_ids_selected is mapped to the nth spot in assigned_locations,
            # the orders of assigned_locations and cell_ids_selected should match.
            assigned_locations = coordinates_data.iloc[index_st_list[process_idx]].iloc[mapped_st_index]\
                                    if index_st_list is not None \
                                    else coordinates_data.iloc[mapped_st_index]
            assigned_locations_list.append(assigned_locations)
            cell_ids_selected = scRNA_norm_data.columns.values[index_sc_list[process_idx]]
            cell_ids_selected_list.append(cell_ids_selected)
        
        cell_ids_selected = np.concatenate(cell_ids_selected_list, axis=0)
        assigned_locations = pd.concat(assigned_locations_list)

    return assigned_locations, cell_ids_selected


def main_cytospace(scRNA_path, cell_type_path,
                   n_cells_per_spot_path, st_cell_type_path, cell_type_fraction_estimation_path=None,
                   spaceranger_path=None, st_path=None, coordinates_path=None,
                   output_folder="cytospace_results", output_prefix="", 
                   mean_cell_numbers=5, downsample_off=False, scRNA_max_transcripts_per_cell=1500,
                   solver_method="lapjv", distance_metric="Pearson_correlation", sampling_method="duplicates",
                   single_cell=False, number_of_selected_spots=10000,
                   sampling_sub_spots=False, number_of_selected_sub_spots=10000,
                   number_of_processors=1, seed=1,
                   plot_off=False, geometry="honeycomb", max_num_cells_plot=50000, num_column=3):
    # For timing execution
    start_time = time.perf_counter()

    # Check paths
    output_path = check_paths(output_folder, output_prefix)

    # Record log
    fout_log = os.path.join(output_path, f"{output_prefix}log.txt")
    with open(fout_log, "w") as f:
        f.write("CytoSPACE log file \n\nStart time: "+str(time.asctime( time.localtime(time.time()) ))+"\n")
        f.write("\nINPUT ARGUMENTS\n")
        f.write("scRNA_path: "+str(scRNA_path)+"\n")
        f.write("cell_type_path: "+str(cell_type_path)+"\n")
        f.write("st_path: "+str(st_path)+"\n")
        f.write("coordinates_path: "+str(coordinates_path)+"\n")
        f.write("n_cells_per_spot_path: "+str(n_cells_per_spot_path)+"\n")
        f.write("cell_type_fraction_estimation_path: "+str(cell_type_fraction_estimation_path)+"\n")
        f.write("st_cell_type_path: "+str(st_cell_type_path)+"\n")
        f.write("output_folder: "+str(output_folder)+"\n")

        f.write("mean_cell_numbers: "+str(mean_cell_numbers)+"\n")
        f.write("downsample_off: "+str(downsample_off)+"\n")
        f.write("scRNA_max_transcripts_per_cell: "+str(scRNA_max_transcripts_per_cell)+"\n")
        f.write("plot_off: "+str(plot_off)+"\n")
        f.write("geometry: "+str(geometry)+"\n")
        f.write("output_prefix: "+str(output_prefix)+"\n")
        f.write("seed: "+str(seed)+"\n")
        f.write("solver_method: "+str(solver_method)+"\n")
        f.write("sampling_method: "+str(sampling_method)+"\n")
        f.write("distance_metric: "+str(distance_metric)+"\n")
        f.write("single_cell: "+str(single_cell)+"\n\n")
        f.write("sampling_sub_spots: "+str(sampling_sub_spots)+"\n\n")

            
    if solver_method == "lapjv" or solver_method == "lapjv_compat":
        solver = import_solver(solver_method)
    else:
        solver = None
    
    # Read data
    print("Read and validate data ...")

    t0 = time.perf_counter()
    scRNA_data, cell_type_data, st_data, coordinates_data, cell_type_fractions_data, n_cells_per_spot_data, st_cell_type_data =\
        read_data(scRNA_path, cell_type_path, 
                  cell_type_fraction_estimation_path, n_cells_per_spot_path, st_cell_type_path,
                  output_path, output_prefix, spaceranger_path, st_path, coordinates_path)
    
    # record all spot IDs
    all_spot_ids = st_data.columns
    
    print(f"Time to read and validate data: {round(time.perf_counter() - t0, 2)} seconds")
    with open(fout_log,"a") as f:
        f.write(f"Time to read and validate data: {round(time.perf_counter() - t0, 2)} seconds\n")
 
    # Set seed
    np.random.seed(seed)
    random.seed(seed)

    t0_core = time.perf_counter()

    ### Parse/compute required inputs
    # cell_number_to_node_assignment
    if single_cell:
        cell_number_to_node_assignment = np.ones(st_data.shape[1]).astype(int)
    else:
        if n_cells_per_spot_data is None:
            print('Estimating number of cells in each spot ...')
            cell_number_to_node_assignment = estimate_cell_number_RNA_reads(st_data, mean_cell_numbers)
            print(f"Time to estimate number of cells per spot: {round(time.perf_counter() - t0_core, 2)} seconds")
        else:
            cell_number_to_node_assignment = n_cells_per_spot_data.values[:, 0].astype(int)

    # cell_type_numbers_int
    if single_cell and (st_cell_type_data is not None):
        cell_types = sorted(st_cell_type_data.iloc[:, 0].unique(), key=str.lower)
        cell_counts = st_cell_type_data.iloc[:, 0].value_counts().reindex(cell_types)
        cell_type_numbers_int = pd.DataFrame(cell_counts)
        cell_type_numbers_int.columns = ['Fraction']
    else:
        if cell_type_fractions_data is None:
            raise ValueError("cell_type_fraction_estimation_path must be specified.")
        # Compute number of cells
        number_of_cells = np.sum(cell_number_to_node_assignment)
        cell_type_numbers_int = get_cell_type_fraction(number_of_cells, cell_type_fractions_data)


    ### Sample scRNA_data
    print("Down/up sample of scRNA-seq data according to estimated cell type fractions")
    t0 = time.perf_counter()

    # subset datasets to intersect_genes
    intersect_genes = st_data.index.intersection(scRNA_data.index)
    scRNA_data_sampled = scRNA_data.loc[intersect_genes, :]
    st_data = st_data.loc[intersect_genes, :]
    

    # downsample scRNA_data_sampled to equal transcript counts per cell
    # so that the assignment is not dependent on expression level
    if not downsample_off:
        scRNA_data_sampled = downsample(scRNA_data_sampled, scRNA_max_transcripts_per_cell)

    with open(fout_log, "a") as f:
        f.write("Number of genes used for mapping: "+str(len(intersect_genes))+"\n")
        f.write("Number of spots satisfying input for mapping: "+str(st_data.shape[1])+"\n")
        f.write("Number of cells satisfying input for mapping: "+str(scRNA_data_sampled.shape[1])+"\n")

    # sample scRNA_data based on cell type composition
    # cell count in scRNA_data_sampled will be equal to cell count (not spot count) in ST data
    # additionally, the cells in scRNA_data_sampled will be in the order of cell types in cell_type_numbers_int
    scRNA_data_sampled =\
        sample_single_cells(scRNA_data_sampled, cell_type_data, cell_type_numbers_int, sampling_method, seed)

    print(f"Time to down/up sample scRNA-seq data: {round(time.perf_counter() - t0, 2)} seconds")

    ### Case-specific parts
    if single_cell:
        with open(fout_log, "a") as f:
            f.write("Number of selected spots: "+str(number_of_selected_spots)+"\n\n")
            f.write("Number of processors: "+str(number_of_processors)+"\n\n")
        
        if st_cell_type_data is not None:
            # recommended way: cell types are specified for each single cell
            index_sc_list = partition_indices(np.arange(scRNA_data_sampled.shape[1]),
                                                split_by_category_list=cell_type_numbers_int.values.flatten(),
                                                split_by_interval_int=number_of_selected_spots,
                                                shuffle=False)
            # order ST indices by cell type (so that e.g., B cells will be at the front and Plasma cells will be at the back)
            # this allows for the cell types of ST partitions to match those of the scRNA partitions above
            index_st_celltype_order = []
            for cell_type in cell_type_numbers_int.index:
                index_st_celltype_order.append(np.where(st_cell_type_data.values[:, 0] == cell_type)[0])
            index_st_celltype_order = np.concatenate(index_st_celltype_order, axis=0)
            index_st_list = partition_indices(index_st_celltype_order,
                                                split_by_category_list=cell_type_numbers_int.values.flatten(),
                                                split_by_interval_int=number_of_selected_spots,
                                                shuffle=False)

            assigned_locations, cell_ids_selected =\
                apply_linear_assignment(scRNA_data_sampled, st_data, coordinates_data, cell_number_to_node_assignment,
                                        solver_method, solver, seed, distance_metric, number_of_processors,
                                        index_sc_list, index_st_list=index_st_list)

        else:

            # cell types are not specified; estimating from cell type fraction data
            index_sc_list = partition_indices(np.arange(scRNA_data_sampled.shape[1]),
                                                split_by_interval_int=number_of_selected_spots,
                                                shuffle=True)
            index_st_list = partition_indices(np.arange(st_data.shape[1]),
                                                split_by_interval_int=number_of_selected_spots,
                                                shuffle=True)

            assigned_locations, cell_ids_selected =\
                apply_linear_assignment(scRNA_data_sampled, st_data, coordinates_data, cell_number_to_node_assignment,
                                        solver_method, solver, seed, distance_metric, number_of_processors,
                                        index_sc_list, index_st_list=index_st_list)

    else:

        if sampling_sub_spots:
            with open(fout_log, "a") as f:
                f.write("Number of selected subspots: "+str(number_of_selected_sub_spots)+"\n\n")
                f.write("Number of processors: "+str(number_of_processors)+"\n\n")
            
            if number_of_selected_sub_spots > np.sum(cell_number_to_node_assignment):
                number_of_selected_sub_spots = np.sum(cell_number_to_node_assignment)

            index_sc_list = partition_indices(np.arange(scRNA_data_sampled.shape[1]),
                                                split_by_interval_int=number_of_selected_sub_spots,
                                                shuffle=True)
            # 1D np.array; repeat each node index by the number of cells in that node
            cell_number_to_node_assignment_aggregate = np.repeat(range(len(cell_number_to_node_assignment)), cell_number_to_node_assignment)
            subsample_spot_idxs_list = partition_indices(cell_number_to_node_assignment_aggregate,
                                                            split_by_interval_int=number_of_selected_sub_spots,
                                                            shuffle=True)
            # aggregate each set of indices in this list back to the (node index) - (count) format
            cell_number_to_node_assignment_list = [np.bincount(subsample_spot_idxs, minlength=(len(cell_number_to_node_assignment)))
                                                    for subsample_spot_idxs in subsample_spot_idxs_list]
                
            assigned_locations, cell_ids_selected =\
                apply_linear_assignment(scRNA_data_sampled, st_data, coordinates_data, cell_number_to_node_assignment,
                                        solver_method, solver, seed, distance_metric, number_of_processors,
                                        index_sc_list, subsampled_cell_number_to_node_assignment_list=cell_number_to_node_assignment_list)
            
        else:

            index_sc_list = partition_indices(np.arange(scRNA_data_sampled.shape[1]), shuffle=False)

            assigned_locations, cell_ids_selected =\
                apply_linear_assignment(scRNA_data_sampled, st_data, coordinates_data, cell_number_to_node_assignment,
                                        solver_method, solver, seed, distance_metric, number_of_processors,
                                        index_sc_list)
                
    
       
    print(f"Total time to run CytoSPACE core algorithm: {round(time.perf_counter() - t0_core, 2)} seconds")
    with open(fout_log,"a") as f:
        f.write(f"Time to run CytoSPACE core algorithm: {round(time.perf_counter() - t0_core, 2)} seconds\n")

    ### Save results
    print('Saving results ...')
    
    # identify unmapped spots
    mapped_spots = assigned_locations.index
    unmapped_spots = np.setdiff1d(list(all_spot_ids), list(mapped_spots)).tolist()

    if len(unmapped_spots) > 0:
        unassigned_locations  = coordinates_data.loc[unmapped_spots]
        unassigned_locations.index = unassigned_locations.index.str.replace("SPOT_", "")
        unassigned_locations["Number of cells"] = 0
        unassigned_locations.to_csv(f"{output_path}/{output_prefix}unassigned_locations.csv", index=True)
        print(f"{len(unmapped_spots)} spots had no cells mapped to them. Saved unfiltered version of assigned locations to {output_path}/{output_prefix}unassigned_locations.csv")

    
    
    save_results(output_path, output_prefix, cell_ids_selected, scRNA_data_sampled if sampling_method == "place_holders" else scRNA_data,
                 assigned_locations, cell_type_data, sampling_method, single_cell)

    if not plot_off:
        if single_cell:
            plot_results(output_path, output_prefix, max_num_cells=max_num_cells_plot, single_cell_ST_mode=True)
        else:
            plot_results(output_path, output_prefix, coordinates_data=coordinates_data, geometry=geometry, num_cols=num_column, max_num_cells=max_num_cells_plot)

    print(f"Total execution time: {round(time.perf_counter() - start_time, 2)} seconds")
    with open(fout_log,"a") as f:
        f.write(f"Total execution time: {round(time.perf_counter() - start_time, 2)} seconds")
    
def run_cytospace():
    arguments = argument_parser()
    main_cytospace(**arguments)
    
if __name__ == "__main__":
    arguments = argument_parser()
    main_cytospace(**arguments)
