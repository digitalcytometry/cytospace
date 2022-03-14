import random
import time
import pandas as pd
import numpy as np

from cytospace.common import read_file, normalize_data, check_paths, argument_parser
from cytospace.post_processing import plot_output, save_results
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
    cell_type_fractions = cell_type_fraction_data.values[0]
    cell_type_numbers = number_of_cells * cell_type_fractions
    cell_type_numbers_int = cell_type_numbers.astype(int)
    number_of_cells_estimated = np.sum(cell_type_numbers_int)
    cell_type_numbers_int[0] += number_of_cells - number_of_cells_estimated

    return cell_type_numbers_int


def solve_linear_assignment_problem(scRNA_data, st_data, cell_type_data,
                                    cell_type_numbers_int, method, coordinates,
                                    cell_number_to_node_assignment, solver_method, solver, seed):
    print("Building cost matrix ...")
    if method == 'shortest_augmenting_path':
        t0 = time.perf_counter()
        distance_repeat, location_repeat, cell_ids_selected, new_cell_index =\
            calculate_cost(scRNA_data, st_data, cell_type_data, cell_type_numbers_int,
                           cell_number_to_node_assignment, seed)
        print(f"Time to build cost matrix: {round(time.perf_counter() - t0, 2)} seconds")

        print('Solving linear assignment problem ...')
        np.random.seed(seed)
        cost_scaled = distance_repeat + 1e-16 * np.random.rand(distance_repeat.shape[0],
                                                               distance_repeat.shape[1])
        t0 = time.perf_counter()
        assignment = call_solver(solver, solver_method, cost_scaled)
        print(f"Time to solve LAP-problem: {round(time.perf_counter() - t0, 2)} seconds")
        assigned_nodes = location_repeat[assignment]
        index = np.transpose(assigned_nodes).tolist()
        assigned_locations = coordinates.iloc[index]

    elif method == 'linear_assignment':
        distance_repeat, location_repeat, cell_ids_selected, new_cell_index =\
            calculate_cost(scRNA_data, st_data, cell_type_data, cell_type_numbers_int,
                           cell_number_to_node_assignment, seed, linear=True)

        print('Solving linear assignment problem ...')
        np.random.seed(seed)
        cost_scaled = 10**6 * distance_repeat + 10 * np.random.rand(distance_repeat.shape[0],
                                                                    distance_repeat.shape[1]) + 1
        cost_scaled = np.transpose(cost_scaled)
        cost_scaled_int = cost_scaled.astype(int)
        cost_scaled_int_list = cost_scaled_int.tolist()
        assignment = match_solution(cost_scaled_int_list)
        assigned_nodes = location_repeat[assignment[:, 0].astype(int)]
        index = assigned_nodes.tolist()
        assigned_locations = coordinates.iloc[index]
    else:
        raise ValueError("Method has to be either 'shortest_augmenting_path' or 'linear_assignment'")

    return assigned_locations, cell_ids_selected, new_cell_index, index, assigned_nodes


def main_cytospace(scRNA_path, cell_type_path, st_path, coordinates_path,
                   cell_type_fraction_estimation_path, output_folder="cytospace_results",
                   method="shortest_augmenting_path", rotation_flag=False, plot_off=False,
                   mean_cell_numbers=5, num_row=3, num_column=4, rotation_degrees=270,
                   output_prefix="", seed=1, delimiter=",", solver_method="lapjv"):
    # For timing execution
    start_time = time.perf_counter()

    # Import LAP-solver based on the solver method
    solver = import_solver(solver_method)

    # Read data
    print("Read and validate data ...")
    t0 = time.perf_counter()
    scRNA_data, cell_type_data, st_data, coordinates_data, cell_type_factions_data =\
        read_data(scRNA_path, cell_type_path, st_path, coordinates_path,
                  cell_type_fraction_estimation_path, delimiter)
    print(f"Time to read and validate data: {round(time.perf_counter() - t0, 2)} seconds")

    # Check paths
    output_path = check_paths(output_folder, output_prefix)

    # Set seed
    random.seed(seed)
    np.random.seed(seed)

    print('Estimating number of cells in each spot ...')
    cell_number_to_node_assignment = estimate_cell_number_RNA_reads(st_data, mean_cell_numbers)

    print('Get cell type fractions ...')
    number_of_cells = np.sum(cell_number_to_node_assignment)
    cell_type_numbers_int = get_cell_type_fraction(number_of_cells, cell_type_factions_data)

    assigned_locations, cell_ids_selected, new_cell_index, index, assigned_nodes =\
        solve_linear_assignment_problem(scRNA_data, st_data, cell_type_data,
                                        cell_type_numbers_int, method, coordinates_data,
                                        cell_number_to_node_assignment, solver_method, solver, seed)

    print('Saving results ...')
    save_results(output_path, output_prefix, cell_ids_selected, assigned_locations,
                 new_cell_index, index, assigned_nodes)

    if not plot_off:
        plot_output(cell_type_fraction_estimation_path, num_row, num_column, rotation_degrees, rotation_flag,
                    output_path, output_prefix, assigned_nodes, new_cell_index, coordinates_data)

    print(f"Total execution time: {round(time.perf_counter() - start_time, 2)} seconds")


def run_cytospace():
    arguments = argument_parser()
    main_cytospace(**arguments)


if __name__ == "__main__":
    arguments = argument_parser()
    main_cytospace(**arguments)
