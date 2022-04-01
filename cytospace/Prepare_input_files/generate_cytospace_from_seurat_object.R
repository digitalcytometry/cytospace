##########################################################
# README
##########################################################
# Script to generate CytoSPACE scRNA-seq and ST input data files  
# given the correspondig Seurat objects.
#
# function inptuts
#
# 1. scRNA Seurat object 
# 2. ST Seurat object 
# 3. Path to the output directory to store the results, i.e., generated CytoSPACE input files
#
##########################################################
library(Seurat)


generate_cytospace_from_seurat_object <- function(scrna_seurat,
                                                   st_seurat,
                                                   fn_out){
  scrna_count <- as.data.frame(as.matrix(GetAssayData(object = scrna_seurat, slot = "counts")))
  scrna_count <- cbind(rownames(scrna_count), scrna_count)
  colnames(scrna_count)[1] <- 'GENES'
  ST_expressions <- as.matrix(st_seurat@assays$Spatial@counts)
  ST_expressions <- cbind(rownames(ST_expressions), ST_expressions)
  colnames(ST_expressions)[1] <- 'GENES'
  coordinates_raw <- st_seurat@images$slice1@coordinates
  coordinates <- coordinates_raw[,-1]
  coordinates <- cbind(rownames(coordinates), coordinates)
  colnames(coordinates)[1] <- 'Spot ID'
  
  print("Writing output to file")
  write.csv(scrna_count, file = paste(fn_out, 'scRNA_data.csv'), row.names = F)
  write.csv(coordinates, file = paste(fn_out, 'Coordinates.csv'), row.names = F)
  write.csv(ST_expressions, file = paste(fn_out, 'ST_data.csv'), row.names = F)
  print("Done")
}

########################### Main #########################
generate_cytospace_from_seurat_object(scrna_seurat,
                                      st_seurat,
                                      fn_out)


