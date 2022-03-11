import pandas as pd
import numpy as np
import scipy
import random

from cytospace.common import normalize_data


def estimate_cell_type_fractions_correlation_based(expressions_st_data, signature_matrix_path,
                                                   correlation_coefficient_limit=0.05):
    # Read data
    expressions_signature_data = pd.read_csv(signature_matrix_path, index_col=0)

    # Find the intersection of genes
    cell_type_name_ordered, index_ordered = np.unique(expressions_signature_data.columns,
                                                      return_index=True)
    expressions_signature_data = expressions_signature_data.iloc[:, index_ordered]
    intersect_genes = expressions_st_data.index.intersection(expressions_signature_data.index)
    expressions_signature_data_intersect_genes = expressions_signature_data.loc[intersect_genes ,:]
    expressions_st_data_intersect_genes = expressions_st_data.loc[intersect_genes ,:]
    expressions_signature = expressions_signature_data_intersect_genes.values.astype(float)
    expressions_ST = expressions_st_data_intersect_genes.values.astype(float)

    # Normalize data
    expressions_tpm_st_log = normalize_data(expressions_ST)
    expressions_tpm_signature_log = normalize_data(expressions_signature)

    # Calculate correlation between ST data and Signature matrix
    no_spots = expressions_tpm_st_log.shape[1]
    no_classes = expressions_tpm_signature_log.shape[1]
    p = np.zeros((no_classes, no_spots))
    r = np.zeros((no_classes, no_spots))
    for i in range(no_classes):
        for j in range(no_spots):
            r[i, j], p[i, j] = scipy.stats.pearsonr(expressions_tpm_signature_log[:, i],
                                                    expressions_tpm_st_log[:, j])

    # Filter out non-significant correlation coefficients
    sum_corr = np.zeros(no_classes)
    for k in range(no_classes):
        index = p[k, :] < correlation_coefficient_limit
        r_cell_type = r[k, :]
        r_cell_type = np.nan_to_num(r_cell_type)
        r_cell_type_significant = r_cell_type[index]
        r_cell_type_significant = r_cell_type_significant - np.min(r_cell_type_significant)
        sum_corr[k] = np.sum(r_cell_type_significant)

    # Estimate cell type fractions
    sum_corr_normalized = sum_corr / np.sum(sum_corr)

    return sum_corr_normalized


def estimate_cell_type_fractions_metagene_based(st_path, marker_genes, seed):
    # Read data and find the intersection of genes
    cell_type_name_ordered, index_ordered = np.unique(marker_genes.columns, return_index=True)
    marker_genes = marker_genes.iloc[:, index_ordered]
    expressions_st_data = pd.read_csv(st_path, header=0, index_col=0)

    expressions_ST = expressions_st_data.values.astype(float)
    expressions_ST = np.nan_to_num(expressions_ST)
    expressions_tpm_ST = (10**6) * (expressions_ST / np.sum(expressions_ST, axis=0, dtype=float))
    expressions_tpm_st_log = np.log2(expressions_tpm_ST + 1)
    expressions_tpm_st_log = np.nan_to_num(expressions_tpm_st_log)
    expressions_tpm_st_log_df = pd.DataFrame(expressions_tpm_st_log,
                                             index=expressions_st_data.index.values,
                                             columns=expressions_st_data.columns.values)

    random.seed(seed)
    noClasses = marker_genes.shape[1]
    noSpots = expressions_st_data.shape[1]
    noGenes = expressions_st_data.shape[0]
    metagene_expression = np.zeros(noClasses)
    for k in range(noClasses):
        # Calculate expression of the marker genes in ST data
        intersect_genes = expressions_tpm_st_log_df.index.intersection(marker_genes.values[: ,k])
        expressions_tpm_st_log_cell_type = expressions_tpm_st_log_df.loc[intersect_genes]
        marker_gene_expression = expressions_tpm_st_log_cell_type.values
        mean_marker_gene_expression = np.mean(marker_gene_expression, axis = 0)

        # Initialize a matrix containing mean expression of random genes
        iteration = 100
        random_expressions = np.zeros((iteration, noSpots))
        for i in range(iteration):
            selected_genes = random.sample(range(noGenes), len(intersect_genes))
            random_expressions[i ,:] = np.mean(expressions_tpm_st_log[selected_genes ,:], axis = 0)

        # Filter out noise based on random chance
        expression_filtered = np.zeros(noSpots)
        for j in range(noSpots):
            t, p = scipy.stats.ttest_ind(mean_marker_gene_expression, random_expressions[: ,j])
            if p < 0.05 and t > 0:
                expression_filtered[j] = mean_marker_gene_expression[j]

        # Estimate cell type fractions
        metagene_expression[k] = np.sum(expression_filtered)
        metagene_expression_normalized = metagene_expression / np.sum(metagene_expression)

    return metagene_expression_normalized