library(Seurat)

##########################################################
# README
##########################################################
# Script to generate CytoSPACE input data files  
# given the correspondig scRNA Seurat object.
#
# function inptuts
#
# 1. scRNA Seurat object 
# 2. Path to the output directory to store the results, i.e., generated CytoSPACE input files
#
##########################################################

generate_cytospace_from_scRNA_seurat_object <- function(scrna_seurat,
                                                        dir_out){
  scrna_count <- as.data.frame(as.matrix(GetAssayData(object = scrna_seurat, slot = "counts")))
  scrna_count_combined <- cbind(rownames(scrna_count), scrna_count)
  colnames(scrna_count_combined)[1] <- 'GENES'
  cell_type_labels <- data.frame(pbmc@meta.data$orig.ident)
  rownames(cell_type_labels) <- colnames(scrna_count)
  cell_type_labels_combined <- cbind(rownames(cell_type_labels), cell_type_labels)
  colnames(cell_type_labels_combined) <- c('Cell IDs', 'CellType')

  print("Writing output to file")
  dir.create(fn_out, showWarnings = FALSE)
  write.csv(scrna_count_combined, file = paste(dir_out, '/scRNA_data.csv', sep = ""), row.names = F)
  write.csv(cell_type_labels_combined, file = paste(dir_out, '/cell_type_labels.csv', sep = ""), row.names = F)
  print("Done")
}

##########################################################
# README
##########################################################
# Script to generate CytoSPACE input data files  
# given the correspondig ST Seurat object.
#
# function inptuts
#
# 1. ST Seurat object 
# 2. Path to the output directory to store the results, i.e., generated CytoSPACE input files
#
##########################################################

generate_cytospace_from_ST_seurat_object <- function(st_seurat,
                                                     dir_out){
  ST_expressions <- as.matrix(st_seurat@assays$Spatial@counts)
  ST_expressions <- cbind(rownames(ST_expressions), ST_expressions)
  colnames(ST_expressions)[1] <- 'GENES'
  coordinates_raw <- st_seurat@images$slice1@coordinates
  coordinates <- coordinates_raw[,-1]
  coordinates <- cbind(rownames(coordinates), coordinates)
  colnames(coordinates)[1] <- 'Spot ID'
  
  print("Writing output to file")
  dir.create(fn_out, showWarnings = FALSE)
  write.csv(coordinates, file = paste(dir_out, '/Coordinates.csv', sep = ""), row.names = F)
  write.csv(ST_expressions, file = paste(dir_out, '/ST_data.csv', sep = ""), row.names = F)
  print("Done")
}
