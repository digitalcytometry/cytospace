##########################################################
# README
##########################################################
# Script to generate CytoSPACE scRNA-seq input data files  
# given the output files of Cell Ranger pipeline.
#
# Input to provide from command line
#
# 1. Path to the directory that contains scRNA data (output of the Cell Ranger pipeline)
# 2. Path to the output directory to store the results, i.e., generated CytoSPACE input files
#
# Run the script from command line with the specified inputs, i.e:
# Rscript generate_cytospace_input_from_cellranger_output.R path_to_scrna_dir_cellranger path_to_output_dir
#
##########################################################
library(hdf5r)
library(rhdf5)
library(Seurat)
library(Matrix)

# Read inputs from command line
args = commandArgs(T)

fin_scrna = args[1]
fn_out = args[2]

barcodes <- as.character(h5read(paste(fin_scrna, 'filtered_feature_bc_matrix.h5', sep =""), "matrix/barcodes"))
h5 <- h5read(paste(fin_scrna, 'filtered_feature_bc_matrix.h5', sep =""), "matrix")
counts <- as.data.frame(as.matrix(sparseMatrix(
  dims = h5$shape,
  i = as.numeric(h5$indices),
  p = as.numeric(h5$indptr),
  x = as.numeric(h5$data),
  index1 = FALSE
)))
colnames(counts) <- barcodes
genes <- as.data.frame(h5[["features"]])$name
non_duplicated_genes <- genes[!duplicated(genes)]
counts_non_duplicated_genes <- counts[!duplicated(genes), ]
rownames(counts_non_duplicated_genes) <- non_duplicated_genes
counts_non_duplicated_genes <- cbind(rownames(counts_non_duplicated_genes), counts_non_duplicated_genes)
colnames(counts_non_duplicated_genes)[1] <- 'GENES'

print("Writing output to file")
write.csv(counts_non_duplicated_genes, file = paste(fn_out, 'scRNA_data.csv'), row.names = F)
print("Done")

sessionInfo()
