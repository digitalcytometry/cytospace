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
# 2. Path to the output directory to store the results. 
#    Default is working directory.
# 3. String to prefix output files with. Default is none.
##########################################################

generate_cytospace_from_scRNA_seurat_object <- function(scrna_seurat,
                                                        dir_out='',fout_prefix='', rna_assay='RNA'){
  scrna_count <- as.data.frame(as.matrix(GetAssayData(object = scrna_seurat, slot = "counts", assay = rna_assay)))
  cell_names <- colnames(scrna_count)
  scrna_count <- cbind(rownames(scrna_count), scrna_count)
  colnames(scrna_count)[1] <- 'GENES'
  cell_type_labels <- data.frame(Idents(scrna_seurat))
  rownames(cell_type_labels) <- cell_names
  cell_type_labels <- cbind(rownames(cell_type_labels), cell_type_labels)
  colnames(cell_type_labels) <- c('Cell IDs', 'CellType')

  print("Writing output to file")
  if(nchar(dir_out)>0){
    dir.create(dir_out, showWarnings = FALSE)
    fout_scrna <- paste0(dir_out,'/',fout_prefix,'scRNA_data.txt')
    fout_labels <- paste0(dir_out,'/',fout_prefix,'cell_type_labels.txt')  
  } else{
    fout_scrna <- paste0(fout_prefix,'scRNA_data.txt')
    fout_labels <- paste0(fout_prefix,'cell_type_labels.txt')  
  }
  write.table(scrna_count, fout_scrna, row.names = F, sep='\t', quote = F)
  write.table(cell_type_labels, file = fout_labels, row.names = F, sep='\t', quote = F)
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
# 2. Path to the output directory to store the results.
#    Default is working directory.
# 3. String to prefix output files with. Default is none.
# 4. Slice name. Default "slice1".
#  
##########################################################

generate_cytospace_from_ST_seurat_object <- function(st_seurat,
                                                     dir_out='',
                                                     fout_prefix='',
                                                     slice='slice1'){
  ST_expressions <- as.matrix(st_seurat@assays$Spatial@counts)
  ST_expressions <- cbind(rownames(ST_expressions), ST_expressions)
  colnames(ST_expressions)[1] <- 'GENES'
  coordinates_raw <- st_seurat@images[[slice]]@coordinates
  coordinates <- coordinates_raw[,c('row','col')]
  coordinates <- cbind(rownames(coordinates), coordinates)
  colnames(coordinates)[1] <- 'Spot ID'
  
  print("Writing output to file")
  if(nchar(dir_out)>0){
    dir.create(dir_out, showWarnings = FALSE)
    fout_st <- paste0(dir_out,'/',fout_prefix,'ST_data.txt')
    fout_coords <- paste0(dir_out,'/',fout_prefix,'Coordinates.txt')  
  } else{
    fout_st <- paste0(fout_prefix,'ST_data.txt')
    fout_coords <- paste0(fout_prefix,'Coordinates.txt')  
  }
  write.table(ST_expressions, fout_st, row.names = F, sep='\t', quote = F)
  write.table(coordinates, file = fout_coords, row.names = F, sep='\t', quote = F)
  print("Done")
}
