import random
import time
import pandas as pd
import numpy as np

from cytospace.common import read_file, normalize_data, check_paths, argument_parser
from cytospace.post_processing import save_results, plot_results
from cytospace.linear_assignment_solvers import (calculate_cost, match_solution, import_solver,
                                                 call_solver)


def read_data(scRNA_path, cell_type_path, st_path, coordinates_path,
              cell_type_fraction_estimation_path, delimiter):
    # Read data
    st_data = read_file(st_path)
    coordinates_data = read_file(coordinates_path)
    scRNA_data = read_file(scRNA_path)
    cell_type_data = read_file(cell_type_path)
    cell_type_faction_data = read_file(cell_type_fraction_estimation_path)

    # Order data to match
    try:
        st_data = st_data[coordinates_data.index]
        scRNA_data = scRNA_data[cell_type_data.index]
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

    return scRNA_data, cell_type_data, st_data, coordinates_data, cell_type_faction_data


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


def solve_linear_assignment_problem(scRNA_data, st_data, cell_type_data,
                                    cell_type_numbers_int, coordinates,
                                    cell_number_to_node_assignment, solver_method, solver, seed):
    distance_repeat, location_repeat, cell_ids_selected, new_cell_index =\
        calculate_cost(scRNA_data, st_data, cell_type_data, cell_type_numbers_int,
                       cell_number_to_node_assignment, seed, solver_method)

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

    return assigned_locations, cell_ids_selected, new_cell_index, index, assigned_nodes



def main_cytospace(scRNA_path, cell_type_path, st_path, coordinates_path,
                   cell_type_fraction_estimation_path, cell_number_estimation_path = None, output_folder="cytospace_results",
                   rotation_flag=True, plot_nonvisium=False, plot_off=False, spot_size=175, plot_marker = 'h',
                   mean_cell_numbers=5, num_row=4, num_column=4, rotation_degrees=270,
                   output_prefix="", seed=1, delimiter=",", solver_method="lapjv"):
    # For timing execution
    start_time = time.perf_counter()

     # Check paths
    output_path = check_paths(output_folder, output_prefix)

    # Record log
    fout_log = output_path / f"{output_prefix}log.txt"
    with open(fout_log, "w") as f:
        f.write("CytoSPACE log file \n\nStart time: "+str(time.asctime( time.localtime(time.time()) ))+"\n")
        f.write("\nINPUT ARGUMENTS\n")
        f.write("scRNA_path: "+str(scRNA_path)+"\n")
        f.write("cell_type_path: "+str(cell_type_path)+"\n")
        f.write("st_path: "+str(st_path)+"\n")
        f.write("coordinates_path: "+str(coordinates_path)+"\n")
        f.write("cell_type_fraction_estimation_path: "+str(cell_type_fraction_estimation_path)+"\n")
        f.write("cell_number_estimation_path: "+str(cell_number_estimation_path)+"\n")
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
        f.write("solver_method: "+str(solver_method)+"\n\n")

    if solver_method == "lapjv" or solver_method == "lapjv_compat":
        solver = import_solver(solver_method)
    else:
        solver = None

    # Read data
    print("Read and validate data ...")
    t0 = time.perf_counter()
    scRNA_data, cell_type_data, st_data, coordinates_data, cell_type_factions_data =\
        read_data(scRNA_path, cell_type_path, st_path, coordinates_path,
                  cell_type_fraction_estimation_path, delimiter)
    print(f"Time to read and validate data: {round(time.perf_counter() - t0, 2)} seconds")
    with open(fout_log,"a") as f:
        f.write(f"Time to read and validate data: {round(time.perf_counter() - t0, 2)} seconds\n")

    # Set seed
    random.seed(seed)
    np.random.seed(seed)

    t0_core = time.perf_counter()
    print('Estimating number of cells in each spot ...')
    if cell_number_estimation_path == None:
         cell_number_to_node_assignment = estimate_cell_number_RNA_reads(st_data, mean_cell_numbers)
    else:
         cell_number_data = pd.read_csv(cell_number_estimation_path, header=0, index_col=0, delimiter=delimiter)
         cell_number_to_node_assignment = cell_number_data.values.astype(int).flatten()
    print(f"Time to estimate number of cells per spot: {round(time.perf_counter() - t0_core, 2)} seconds")

    print('Get cell type fractions ...')
    number_of_cells = np.sum(cell_number_to_node_assignment)
    cell_type_numbers_int = get_cell_type_fraction(number_of_cells, cell_type_factions_data)

    assigned_locations, cell_ids_selected, new_cell_index, index, assigned_nodes =\
        solve_linear_assignment_problem(scRNA_data, st_data, cell_type_data,
                                        cell_type_numbers_int, coordinates_data,
                                        cell_number_to_node_assignment, solver_method, solver, seed)
    print(f"Total time to run CytoSPACE core algorithm: {round(time.perf_counter() - t0_core, 2)} seconds")
    with open(fout_log,"a") as f:
        f.write(f"Time to run CytoSPACE core algorithm: {round(time.perf_counter() - t0_core, 2)} seconds\n")

    print('Saving results ...')
    assigned_locations_path = str(output_path / f'{output_prefix}assigned_locations.csv')
    save_results(output_path, output_prefix, cell_ids_selected, assigned_locations,
                 new_cell_index, index, assigned_nodes, st_path, coordinates_path,
                 cell_type_path, assigned_locations_path)

    if not plot_off:
        output_filename = output_path / "plot_cell_type_locations.pdf"
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
