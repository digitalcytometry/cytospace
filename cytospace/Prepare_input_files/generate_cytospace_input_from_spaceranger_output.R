##########################################################
# README
##########################################################
# Script to generate ST input data files  
# given the output files of the Space Ranger (10x Visium) pipeline.
#
# Input to provide from command line
#
# 1. Path to the directory that contains ST data (output of the Space Ranger pipeline)
# 2. Path to the output directory to store the results, i.e., generated CytoSPACE input files
#
# Run the script from command line with the specified inputs, i.e:
# Rscript generate_cytospace_input_from_spaceranger_output.R path_to_ST_dir_spaceranger path_to_output_dir
# 
# If writing the output as sparse matrix instead of CSV file, you can append --sparse to the end of the command:
# Rscript generate_cytospace_input_from_spaceranger_output.R path_to_ST_dir_spaceranger path_to_output_dir --sparse
##########################################################
library(hdf5r)
library(Seurat)
library(Matrix)

# Read inputs from command line
args = commandArgs(T)

fin_ST = args[1]
fn_out = args[2]
write_sparse = ifelse(length(args) > 2, args[3] == '--sparse', FALSE)

ST_data <- Load10X_Spatial(fin_ST)

ST_expressions <- ST_data@assays$Spatial@counts
spot_names <- colnames(ST_expressions)
gene_names <- rownames(ST_expressions)

if (!dir.exists(fn_out)) {
  dir.create(fn_out, showWarnings = FALSE)
}

if (write_sparse) {
  fout_st <- file.path(fn_out, 'ST_data.mtx')
  fout_genes <- file.path(fn_out, 'ST_data_genes.tsv')
  fout_spots <- file.path(fn_out, 'ST_data_cells.tsv')

  Matrix::writeMM(ST_expressions, fout_st)
  write.table(as.data.frame(gene_names), fout_genes, row.names = F, col.names = F, sep='\t', quote = F)
  write.table(as.data.frame(spot_names), fout_spots, row.names = F, col.names = F, sep='\t', quote = F)
} else {
  ST_expressions <- as.data.frame(as.matrix(ST_expressions))
  ST_expressions <- cbind(rownames(ST_expressions), ST_expressions)
  colnames(ST_expressions)[1] <- 'GENES'

  print("Writing output to file")
  fout_st <- file.path(fn_out, 'ST_data.csv')
  write.csv(ST_expressions, fout_st, row.names = F, quote = F)
}

# write cell type labels file
coordinates <- ST_data@images$slice1@coordinates[, -1]
coordinates <- cbind(rownames(coordinates), coordinates)
colnames(coordinates)[1] <- 'Spot ID'

fout_coords <- file.path(fn_out, 'Coordinates.csv')
write.csv(coordinates, fout_coords, row.names = F, quote = F)

print("Done")
