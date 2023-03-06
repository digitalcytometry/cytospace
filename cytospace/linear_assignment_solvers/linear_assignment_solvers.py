import numpy as np
import pandas as pd
import random
import time
from ortools.graph import pywrapgraph
from cytospace.common import normalize_data, matrix_correlation_pearson, matrix_correlation_spearman
from scipy.spatial import distance


def import_solver(solver_method):
    try:
        if solver_method == "lapjv_compat":
            from lap import lapjv
            solver = lapjv
        elif solver_method == "lapjv":
            from lapjv import lapjv
            solver = lapjv
        else:
            raise NotImplementedError(f"The solver {solver_method} is not a supported solver "
                                      "for the shortest augmenting path method, choose between "
                                      "'lapjv' and 'lapjv_compat'.")
    except ModuleNotFoundError:
        raise ModuleNotFoundError("The Python package containing the solver_method option "
                                  f"you have chosen {solver_method} was not found. If you "
                                  "selected 'lapjv_compat' solver, install package 'lap'"
                                  "by running 'pip install lap==0.4.0'. If you selected 'lapjv'"
                                  "solver, install package 'lapjv' by running `pip intall lapjv==1.3.14'"
                                  "or check the package home page for further instructions.")

    return solver


def call_solver(solver, solver_method, cost_scaled):
    if solver_method == "lapjv_compat":
        _, _, y = solver(cost_scaled)
    elif solver_method == "lapjv":
        _, y, _ = solver(cost_scaled)

    return y

def calculate_cost(expressions_tpm_scRNA_log, expressions_tpm_st_log, cell_number_to_node_assignment,
                    solver_method, distance_metric):
    print("Building cost matrix ...")
    t0 = time.perf_counter()
    if solver_method=="lap_CSPR":
        if distance_metric=="Pearson_correlation":
           cost = -np.transpose(matrix_correlation_pearson(expressions_tpm_st_log, expressions_tpm_scRNA_log))
        elif distance_metric=="Spearman_correlation":
           cost = -np.transpose(matrix_correlation_spearman(expressions_tpm_st_log, expressions_tpm_scRNA_log))
        elif distance_metric=="Euclidean":
           cost = np.transpose(distance.cdist(np.transpose(expressions_tpm_scRNA_log), np.transpose(expressions_tpm_st_log), 'euclidean'))
    else:
        if distance_metric=="Pearson_correlation":
           cost = -matrix_correlation_pearson(expressions_tpm_scRNA_log, expressions_tpm_st_log)
        elif distance_metric=="Spearman_correlation":
           cost = -matrix_correlation_spearman(expressions_tpm_scRNA_log, expressions_tpm_st_log)
        elif distance_metric=="Euclidean":
           cost = np.transpose(distance.cdist(np.transpose(expressions_tpm_scRNA_log), np.transpose(expressions_tpm_st_log), 'euclidean'))

    location_repeat = np.zeros(cost.shape[1])
    counter = 0
    location_repeat = np.repeat(np.arange(len(cell_number_to_node_assignment)), cell_number_to_node_assignment)

    location_repeat = location_repeat.astype(int)
    distance_repeat = cost[location_repeat, :]
    print(f"Time to build cost matrix: {round(time.perf_counter() - t0, 2)} seconds")

    return distance_repeat, location_repeat


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
