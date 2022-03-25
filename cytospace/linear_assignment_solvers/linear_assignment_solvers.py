import numpy as np
import random
from ortools.graph import pywrapgraph

from cytospace.common import normalize_data, matrix_correlation


def import_solver(solver_method):
    try:
        if solver_method == "lap":
            from lap import lapjv
            solver = lapjv
        elif solver_method == "lapjv":
            from lapjv import lapjv
            solver = lapjv
        elif solver_method == "lapsolver":
            from lapsolver import solve_dense
            solver = solve_dense
        else:
            raise NotImplementedError(f"The method {solver_method} is not an supported method, "
                                      "chose between lap, lapjv, and lapsolver.")
    except ModuleNotFoundError:
        raise ModuleNotFoundError(f"The solver_method option you have chosen {solver_method} "
                                  "depends on the python package with the same name. Please "
                                  f"install it by running ´pip install {solver_method}´ or confer "
                                  "with the package home page for further instructions.")

    return solver


def call_solver(solver, solver_method, cost_scaled):
    if solver_method == "lap":
        _, _, y = solver(cost_scaled)
    elif solver_method == "lapjv":
        _, y, _ = solver(cost_scaled)
    elif solver_method == "lapsolver":
        _, y = solver(cost_scaled)

    return y


def calculate_cost(expressions_scRNA_data, expressions_st_data, cell_type_labels, cell_type_numbers_int,
                   cell_number_to_node_assignment, seed, linear=False):
    # Find intersection genes
    intersect_genes = expressions_st_data.index.intersection(expressions_scRNA_data.index)
    expressions_scRNA_data_intersect_genes = expressions_scRNA_data.loc[intersect_genes, :]
    expressions_st_data_intersect_genes = expressions_st_data.loc[intersect_genes, :]
    expressions_scRNA = expressions_scRNA_data_intersect_genes.values.astype(float)
    expressions_st = expressions_st_data_intersect_genes.values.astype(float)

    # Data normalization
    expressions_tpm_st_log = normalize_data(expressions_st)
    expressions_tpm_scRNA_log = normalize_data(expressions_scRNA)

    # Down/up sample of scRNA-seq data according to estimated cell type fractions
    unique_cell_type_labels = sorted(np.unique(cell_type_labels.values), key=str.lower)
    number_classes = len(unique_cell_type_labels)
    new_cell_type = np.zeros(number_classes + 1)

    # Build cost matrix
    np.random.seed(seed)
    random.seed(seed)
    sampled_index_total = []
    for k in range(0, number_classes):
        cell_type_index = np.nonzero(cell_type_labels.values == unique_cell_type_labels[k])[0].tolist()
        fractions_cells = len(cell_type_index)
        fractions_beads = int(cell_type_numbers_int[k])
        if fractions_beads > fractions_cells:
            sampled_index = np.random.choice(cell_type_index,
                                             fractions_beads - fractions_cells).tolist()
            new_cells = np.concatenate((expressions_tpm_scRNA_log[:, cell_type_index],
                                        expressions_tpm_scRNA_log[:, sampled_index]), axis=1)
            sampled_index_total += cell_type_index + sampled_index
        else:
            sampled_index = random.sample(cell_type_index, fractions_beads)
            new_cells = expressions_tpm_scRNA_log[:, sampled_index]
            sampled_index_total += sampled_index

        if k == 0:
            sampled_cells = new_cells
        else:
            sampled_cells = np.concatenate((sampled_cells, new_cells), axis=1)

        new_cell_type[k + 1] = new_cell_type[k] + new_cells.shape[1]

    if linear:
        cost = -np.transpose(matrix_correlation(expressions_tpm_st_log, sampled_cells))
    else:
        cost = -matrix_correlation(sampled_cells, expressions_tpm_st_log)

    location_repeat = np.zeros(cost.shape[1])
    counter = 0
    for value, repeat in enumerate(cell_number_to_node_assignment):
        location_repeat[counter: counter + repeat] = value
        counter += repeat

    location_repeat = location_repeat.astype(int)
    distance_repeat = cost[location_repeat, :]
    cell_ids = expressions_scRNA_data.columns.values
    cell_ids_selected = cell_ids[sampled_index_total]

    return distance_repeat, location_repeat, cell_ids_selected, new_cell_type


def match_solution(cost):
    rows = len(cost)
    cols = len(cost[0])
    assignment_mat = np.zeros((rows, 2))
    assignment = pywrapgraph.LinearSumAssignment()
    for worker in range(rows):
        for task in range(cols):
            if cost[worker][task]:
                assignment.AddArcWithCost(worker, task, cost[worker][task])

    solve_status = assignment.Solve()
    if solve_status == assignment.OPTIMAL:
        print('Total cost = ', assignment.OptimalCost())
        print()
        for i in range(0, assignment.NumNodes()):
            assignment_mat[i, 0] = assignment.RightMate(i)
            assignment_mat[i, 1] = assignment.AssignmentCost(i)
    elif solve_status == assignment.INFEASIBLE:
        print('No assignment is possible.')
    elif solve_status == assignment.POSSIBLE_OVERFLOW:
        print('Some input costs are too large and may cause an integer overflow.')
    else:
        raise ValueError("The assignment failed")

    return assignment_mat
