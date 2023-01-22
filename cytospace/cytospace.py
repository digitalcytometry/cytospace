import random
import time
import pandas as pd
import numpy as np
import concurrent.futures
from pandas.core.common import flatten

from cytospace.common import read_file, normalize_data, check_paths, argument_parser
from cytospace.post_processing import save_results, plot_results
from cytospace.linear_assignment_solvers import (calculate_cost, match_solution, import_solver,
                                                 call_solver)


def read_data(scRNA_path, cell_type_path, st_path, coordinates_path,
              cell_type_fraction_estimation_path, n_cells_per_spot_path, delimiter):
    # Read data
    st_data = read_file(st_path)
    coordinates_data = read_file(coordinates_path)
    scRNA_data = read_file(scRNA_path)
    cell_type_data = read_file(cell_type_path)
    cell_type_faction_data = read_file(cell_type_fraction_estimation_path)
    if n_cells_per_spot_path is not None:
        n_cells_per_spot_data = read_file(n_cells_per_spot_path)
    else:
        n_cells_per_spot_data = None

    # Order data to match
    try:
        st_data = st_data[coordinates_data.index]
        scRNA_data = scRNA_data[cell_type_data.index]
        if n_cells_per_spot_data is not None:
            n_cells_per_spot_data = n_cells_per_spot_data.transpose(copy=False)
            n_cells_per_spot_data = n_cells_per_spot_data[coordinates_data.index]
            n_cells_per_spot_data = n_cells_per_spot_data.transpose(copy=False)
    except Exception as e:
        raise IndexError("The ST data: {st_path} and coordinates data: {coordinates_path} have to "
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

    return scRNA_data, cell_type_data, st_data, coordinates_data, cell_type_faction_data, n_cells_per_spot_data


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
    cell_type_fraction_data = cell_type_fraction_data.reindex(sorted(cell_type_fraction_data.columns, key=str.lower), axis = 1)
    cell_type_fractions = cell_type_fraction_data.values[0]
    cell_type_numbers = number_of_cells * cell_type_fractions
    cell_type_numbers_int = cell_type_numbers.astype(int)
    number_of_cells_estimated = np.sum(cell_type_numbers_int)
    cell_type_numbers_int[0] += number_of_cells - number_of_cells_estimated

    return cell_type_numbers_int


def solve_linear_assignment_problem(scRNA_data, st_data, cell_type_data, cell_type_numbers_int, coordinates,
                                    cell_number_to_node_assignment, solver_method, sampling_method, solver,
                                    seed, distance_metric):
    distance_repeat, location_repeat, cell_ids_selected, new_cell_index, cell_ids_new, all_cells_save =\
        calculate_cost(scRNA_data, st_data, cell_type_data, cell_type_numbers_int,
                       cell_number_to_node_assignment, seed, solver_method, sampling_method, distance_metric)

    if solver_method == 'lapjv' or solver_method == 'lapjv_compat':
        print('Solving linear assignment problem ...')
        np.random.seed(seed)
        cost_scaled = distance_repeat + 1e-16 * np.random.rand(distance_repeat.shape[0],
                                                               distance_repeat.shape[1])
        t0 = time.perf_counter()
        assignment = call_solver(solver, solver_method, cost_scaled)
        print(f"Time to solve linear assignment problem: {round(time.perf_counter() - t0, 2)} seconds")
        assigned_nodes = location_repeat[assignment]
        index = np.transpose(assigned_nodes).tolist()
        assigned_locations = coordinates.iloc[index]

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
        index = assigned_nodes.tolist()
        assigned_locations = coordinates.iloc[index]
    else:
        raise ValueError("Invalid solver_method provided")

    return assigned_locations, cell_ids_selected, new_cell_index, index, assigned_nodes, cell_ids_new, all_cells_save


def apply_linear_assignment(index_total,cell_type_factions_data,number_of_selected_cells,scRNA_data, st_data, cell_type_data,
                                    coordinates_data, cell_number_to_node_assignment, solver_method, sampling_method, solver, seed, distance_metric):

    expressions_st_selected = st_data.iloc[:, index_total]
    selected_index_sc = list(range(scRNA_data.shape[1]))
    index_sc = np.random.choice(selected_index_sc, number_of_selected_cells).tolist()    
    expressions_sc_selected = scRNA_data.iloc[:, index_sc]
    cell_type_data_selected = cell_type_data.iloc[index_sc]
    cell_number_to_node_assignment_selected = cell_number_to_node_assignment[index_total]
    coordinates_data_selected = coordinates_data.iloc[index_total,:]

    print('Get cell type fractions ...')
    number_of_cells = np.sum(cell_number_to_node_assignment_selected)
    cell_type_numbers_int_selected = get_cell_type_fraction(number_of_cells, cell_type_factions_data)
    assigned_locations, cell_ids_selected, new_cell_index, index, assigned_nodes, cell_ids_new, all_cells_save =\
    solve_linear_assignment_problem(expressions_sc_selected, expressions_st_selected, cell_type_data_selected,
                                    cell_type_numbers_int_selected, coordinates_data_selected,
                                    cell_number_to_node_assignment_selected, solver_method, sampling_method, solver, seed, distance_metric)
    return assigned_locations, cell_ids_selected, new_cell_index, index, assigned_nodes, cell_ids_new, all_cells_save


def main_cytospace(scRNA_path, cell_type_path, st_path, coordinates_path,
                   cell_type_fraction_estimation_path, n_cells_per_spot_path, output_folder="cytospace_results",
                   rotation_flag=True, plot_nonvisium=False, plot_off=False, spot_size=175, plot_marker = 'h',
                   mean_cell_numbers=5, num_row=4, num_column=4, rotation_degrees=270,
                   output_prefix="", seed=1, delimiter=",", solver_method="lapjv", sampling_method="duplicates",
                   distance_metric="Pearson_correlation", number_of_selected_cells=10000, number_of_selected_spots=10000,
                   number_of_processors=4, single_cell=False, sampling_sub_spots=False, number_of_selected_sub_spots=10000):
    # For timing execution
    start_time = time.perf_counter()

    # Check paths
    output_path = check_paths(output_folder, output_prefix)
    assigned_locations_path = str(output_path / f'{output_prefix}assigned_locations.csv')

    # Record log
    fout_log = output_path / f"{output_prefix}log.txt"
    with open(fout_log, "w") as f:
        f.write("CytoSPACE log file \n\nStart time: "+str(time.asctime( time.localtime(time.time()) ))+"\n")
        f.write("\nINPUT ARGUMENTS\n")
        f.write("scRNA_path: "+str(scRNA_path)+"\n")
        f.write("cell_type_path: "+str(cell_type_path)+"\n")
        f.write("st_path: "+str(st_path)+"\n")
        f.write("coordinates_path: "+str(coordinates_path)+"\n")
        f.write("n_cells_per_spot_path: "+str(n_cells_per_spot_path)+"\n")
        f.write("cell_type_fraction_estimation_path: "+str(cell_type_fraction_estimation_path)+"\n")
        f.write("output_folder: "+str(output_folder)+"\n")
        f.write("rotation_flag: "+str(rotation_flag)+"\n")
        f.write("plot_nonvisium: "+str(plot_nonvisium)+"\n")
        f.write("plot_off: "+str(plot_off)+"\n")
        f.write("spot_size: "+str(spot_size)+"\n")
        f.write("plot_marker: "+str(plot_marker)+"\n")
        f.write("mean_cell_numbers: "+str(mean_cell_numbers)+"\n")
        f.write("num_row: "+str(num_row)+"\n")
        f.write("num_column: "+str(num_column)+"\n")
        f.write("rotation_degrees: "+str(rotation_degrees)+"\n")
        f.write("output_prefix: "+str(output_prefix)+"\n")
        f.write("seed: "+str(seed)+"\n")
        f.write("delimiter: "+str(delimiter)+"\n")
        f.write("solver_method: "+str(solver_method)+"\n")
        f.write("sampling_method: "+str(sampling_method)+"\n")
        f.write("distance_metric: "+str(distance_metric)+"\n")
        f.write("single_cell: "+str(single_cell)+"\n\n")
        f.write("sampling_sub_spots: "+str(sampling_sub_spots)+"\n\n")
        f.write("number_of_selected_sub_spots: "+str(number_of_selected_sub_spots)+"\n\n")

            
    if solver_method == "lapjv" or solver_method == "lapjv_compat":
        solver = import_solver(solver_method)
    else:
        solver = None
    
    # Read data
    print("Read and validate data ...")
    t0 = time.perf_counter()
    scRNA_data, cell_type_data, st_data, coordinates_data, cell_type_factions_data, n_cells_per_spot_data =\
        read_data(scRNA_path, cell_type_path, st_path, coordinates_path,
                  cell_type_fraction_estimation_path, n_cells_per_spot_path, delimiter)
    print(f"Time to read and validate data: {round(time.perf_counter() - t0, 2)} seconds")
    with open(fout_log,"a") as f:
        f.write(f"Time to read and validate data: {round(time.perf_counter() - t0, 2)} seconds\n")
 
    # Set seed
    random.seed(seed)
    np.random.seed(seed)

    t0_core = time.perf_counter()

    if single_cell:
        
        size_of_ST_data = st_data.shape[1]

        # Validate "number of selected spots" input
        if size_of_ST_data < number_of_selected_spots:
            number_of_selected_spots = size_of_ST_data
            print("Since number_of_selected_spots is higher than the size of ST data, number_of_selected_spots has been set to the size of ST data")
        
        size_of_scRNA_data = scRNA_data.shape[1]
        
        # Validate "number of selected spots" input
        if size_of_scRNA_data < number_of_selected_cells:
            number_of_selected_cells = size_of_scRNA_data
            print("Since number_of_selected_cells is higher than the size of scRNA-seq data, number_of_selected_cells has been set to the size of scRNA-seq data")
        
        with open(fout_log, "a") as f:
            f.write("number_of_selected_cells: "+str(number_of_selected_cells)+"\n\n")
            f.write("number_of_selected_spots: "+str(number_of_selected_spots)+"\n\n")
            f.write("number_of_processors: "+str(number_of_processors)+"\n\n")
        
        cell_number_to_node_assignment = np.ones(st_data.shape[1]).astype(int)

        max_value = int(size_of_ST_data/(number_of_selected_spots + 1)) + 1
        index_total = [[]*number_of_selected_spots]*max_value
        selected_index_st = list(range(size_of_ST_data))
        all_cells_save_combined = pd.DataFrame()
        cell_ids_selected_combined = pd.DataFrame()
        assigned_locations_combined = pd.DataFrame()
        index_combined = []
    
        print(f"Number of required processors: {max_value}")
    
        for i in range(max_value):
            index_total[i] = np.random.choice(selected_index_st, min(number_of_selected_spots,len(selected_index_st)), replace = False).tolist()
            [selected_index_st.remove(index_total[i][j]) for j in range(len(index_total[i]))]
    
        iter_number = 1     
        if max_value <= number_of_processors:
    
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = [executor.submit(apply_linear_assignment, index_total[i], cell_type_factions_data, number_of_selected_cells, scRNA_data, st_data, cell_type_data, coordinates_data,
                                            cell_number_to_node_assignment, solver_method, sampling_method, solver, seed, distance_metric) for i in range(max_value)]
                for f in concurrent.futures.as_completed(results):
                    assigned_locations, cell_ids_selected, new_cell_index, index, assigned_nodes, cell_ids_new, all_cells_save = f.result()
    
                    all_cells_save_combined = pd.concat([all_cells_save_combined, all_cells_save], axis=1)
                    cell_ids_selected_combined = pd.concat([cell_ids_selected_combined, pd.DataFrame(cell_ids_selected)])
                    assigned_locations_combined = pd.concat([assigned_locations_combined, assigned_locations])
                    index_combined = index_combined + index
                print(f"Iteration: {iter_number}")
                iter_number = iter_number + 1
    
        else:
            counter = 0
            while max_value > 0:
    
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    results = [executor.submit(apply_linear_assignment, index_total[i], cell_type_factions_data, number_of_selected_cells, scRNA_data, st_data, cell_type_data, coordinates_data,
                                                cell_number_to_node_assignment, solver_method, sampling_method, solver, seed, distance_metric) for i in range(counter, counter + min(max_value,number_of_processors))]
                    for f in concurrent.futures.as_completed(results):
                        assigned_locations, cell_ids_selected, new_cell_index, index, assigned_nodes, cell_ids_new, all_cells_save = f.result()
    
                        all_cells_save_combined = pd.concat([all_cells_save_combined, all_cells_save], axis=1)
                        cell_ids_selected_combined = pd.concat([cell_ids_selected_combined, pd.DataFrame(cell_ids_selected)])
                        assigned_locations_combined = pd.concat([assigned_locations_combined, assigned_locations])
                        index_combined = index_combined + index
                        
                    counter = counter + min(max_value,number_of_processors)
                    max_value = max_value - number_of_processors
                    print(f"Iteration: {iter_number}")
                    iter_number = iter_number + 1
    
        all_cells_save = all_cells_save_combined
        cell_ids_new = all_cells_save.columns
        cell_ids_selected = cell_ids_selected_combined
        cell_ids_selected = cell_ids_selected.values.tolist()
        cell_ids_selected = flatten(cell_ids_selected)
        assigned_locations = assigned_locations_combined
        index = index_combined
                                     
    else:
        
        if n_cells_per_spot_data is None:
            print('Estimating number of cells in each spot ...')
            cell_number_to_node_assignment = estimate_cell_number_RNA_reads(st_data, mean_cell_numbers)
            print(f"Time to estimate number of cells per spot: {round(time.perf_counter() - t0_core, 2)} seconds")
        else:
            cell_number_to_node_assignment = n_cells_per_spot_data.values[:, 0].astype(int)


        if  sampling_sub_spots:
            if number_of_selected_sub_spots > np.sum(cell_number_to_node_assignment):
                number_of_selected_sub_spots = np.sum(cell_number_to_node_assignment)
                
            cell_number_to_node_assignment_aggregate = np.zeros((np.sum(cell_number_to_node_assignment),1))
            counter = 0
            for i in range(len(cell_number_to_node_assignment)):
                cell_number_to_node_assignment_aggregate[counter:counter + cell_number_to_node_assignment[i]] = i
                counter = counter + cell_number_to_node_assignment[i]
    
            index_sub_spot = np.random.choice(range(np.sum(cell_number_to_node_assignment)), number_of_selected_sub_spots).tolist()    
            cell_number_to_node_assignment_aggregate_selected = cell_number_to_node_assignment_aggregate[index_sub_spot]
            for i in range(len(cell_number_to_node_assignment)):
                cell_number_to_node_assignment_aggregate_selected_bool = cell_number_to_node_assignment_aggregate_selected == i
                cell_number_to_node_assignment[i] = sum(bool(x) for x in cell_number_to_node_assignment_aggregate_selected_bool)
            
        
        print('Get cell type fractions ...')
        number_of_cells = np.sum(cell_number_to_node_assignment)
        cell_type_numbers_int = get_cell_type_fraction(number_of_cells, cell_type_factions_data)
    
        assigned_locations, cell_ids_selected, new_cell_index, index, assigned_nodes, cell_ids_new, all_cells_save =\
            solve_linear_assignment_problem(scRNA_data, st_data, cell_type_data,
                                            cell_type_numbers_int, coordinates_data,
                                            cell_number_to_node_assignment, solver_method, sampling_method, solver, seed, distance_metric)
            
            
    print(f"Total time to run CytoSPACE core algorithm: {round(time.perf_counter() - t0_core, 2)} seconds")
    with open(fout_log,"a") as f:
        f.write(f"Time to run CytoSPACE core algorithm: {round(time.perf_counter() - t0_core, 2)} seconds\n")

    print('Saving results ...')
    save_results(output_path, output_prefix, cell_ids_selected, cell_ids_new, all_cells_save, assigned_locations,
                 new_cell_index, index, assigned_nodes, st_path, coordinates_path,
                 cell_type_path, assigned_locations_path, single_cell, sampling_method)

    if not plot_off and not single_cell:
        output_filename = output_path / (output_prefix + "plot_cell_type_locations.pdf")
        plot_results(assigned_locations_path, coordinates_path, output_filename, num_row, num_column,
                     rotation_flag, plot_nonvisium, rotation_degrees, spot_size, plot_marker)    

    print(f"Total execution time: {round(time.perf_counter() - start_time, 2)} seconds")
    with open(fout_log,"a") as f:
        f.write(f"Total execution time: {round(time.perf_counter() - start_time, 2)} seconds\n")
    
def run_cytospace():
    arguments = argument_parser()
    main_cytospace(**arguments)
    
if __name__ == "__main__":
    arguments = argument_parser()
    main_cytospace(**arguments)
