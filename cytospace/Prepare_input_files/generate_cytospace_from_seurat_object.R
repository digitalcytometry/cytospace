library(Seurat)
library(Matrix)

##########################################################
# README
##########################################################
# Script to generate CytoSPACE input data files  
# given the corresponding scRNA Seurat object.
#
# function inputs
#
# 1. scRNA Seurat object 
# 2. Path to the output directory to store the results. 
#    Default is working directory.
# 3. String to prefix output files with. Default is none.
# 4. Assay to pull counts from. Default is RNA.
##########################################################

generate_cytospace_from_scRNA_seurat_object <- function(scrna_seurat,
                                                        dir_out='.', fout_prefix='',
                                                        write_sparse=FALSE, rna_assay='RNA'){
  scrna_count <- GetAssayData(scrna_seurat, slot='counts', assay=rna_assay)
  cell_names <- colnames(scrna_count)
  gene_names <- rownames(scrna_count)

  if (dir_out == '') {
    # compatibility with rhe previous version of this script
    message('Note: the output directory will be set to the current working directory.')
  }
  if (!dir.exists(dir_out)) {
    dir.create(dir_out, showWarnings = FALSE)
  }

  if (write_sparse) {
    fout_scrna <- file.path(dir_out, paste0(fout_prefix, 'scRNA_data.mtx'))
    fout_genes <- file.path(dir_out, paste0(fout_prefix, 'genes.tsv'))
    fout_cells <- file.path(dir_out, paste0(fout_prefix, 'cells.tsv'))

    Matrix::writeMM(scrna_count, fout_scrna)
    write.table(as.data.frame(gene_names), fout_genes, row.names = F, col.names = F, sep='\t', quote = F)
    write.table(as.data.frame(cell_names), fout_cells, row.names = F, col.names = F, sep='\t', quote = F)
  } else {
    scrna_count <- as.data.frame(as.matrix(scrna_count))
    scrna_count <- cbind(rownames(scrna_count), scrna_count)
    colnames(scrna_count)[1] <- 'GENES'

    print("Writing output to file")
    fout_scrna <- file.path(dir_out, paste0(fout_prefix, 'scRNA_data.txt'))
    write.table(scrna_count, fout_scrna, row.names = F, sep='\t', quote = F)
  }

  # write cell type labels file
  cell_type_labels <- data.frame(Idents(scrna_seurat))
  rownames(cell_type_labels) <- cell_names
  cell_type_labels <- cbind(rownames(cell_type_labels), cell_type_labels)
  colnames(cell_type_labels) <- c('Cell IDs', 'CellType')

  fout_labels <- file.path(dir_out, paste0(fout_prefix, 'cell_type_labels.txt'))
  write.table(cell_type_labels, fout_labels, row.names = F, sep='\t', quote = F)

  print("Done")
}


##########################################################
# README
##########################################################
# Script to generate CytoSPACE input data files  
# given the corresponding ST Seurat object.
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
                                                     dir_out='.', fout_prefix='',
                                                     write_sparse=FALSE, slice='slice1'){
  ST_expressions <- st_seurat@assays$Spatial@counts
  spot_names <- colnames(ST_expressions)
  gene_names <- rownames(ST_expressions)

  if (dir_out == '') {
    # compatibility with the previous version of this script
    message('Note: the output directory will be set to the current working directory.')
  }
  if (!dir.exists(dir_out)) {
    dir.create(dir_out, showWarnings = FALSE)
  }

  if (write_sparse) {
    fout_st <- file.path(dir_out, paste0(fout_prefix, 'ST_data.mtx'))
    fout_genes <- file.path(dir_out, paste0(fout_prefix, 'genes.tsv'))
    fout_spots <- file.path(dir_out, paste0(fout_prefix, 'cells.tsv'))

    Matrix::writeMM(ST_expressions, fout_st)
    write.table(as.data.frame(gene_names), fout_genes, row.names = F, col.names = F, sep='\t', quote = F)
    write.table(as.data.frame(spot_names), fout_spots, row.names = F, col.names = F, sep='\t', quote = F)
  } else {
    ST_expressions <- as.data.frame(as.matrix(ST_expressions))
    ST_expressions <- cbind(rownames(ST_expressions), ST_expressions)
    colnames(ST_expressions)[1] <- 'GENES'

    print("Writing output to file")
    fout_st <- file.path(dir_out, paste0(fout_prefix, 'ST_data.txt'))
    write.table(ST_expressions, fout_st, row.names = F, sep='\t', quote = F)
  }

  # write cell type labels file
  coordinates <- st_seurat@images[[slice]]@coordinates[, c('row', 'col')]
  coordinates <- cbind(rownames(coordinates), coordinates)
  colnames(coordinates)[1] <- 'Spot ID'

  fout_coords <- file.path(dir_out, paste0(fout_prefix, 'Coordinates.txt'))
  write.table(coordinates, fout_coords, row.names = F, sep='\t', quote = F)

  print("Done")
}