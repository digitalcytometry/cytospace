
##########################################################
# README
##########################################################
# Script to compute overall cell type fractions for ST
# sample given scRNA set with cell type annotations as
# well as fractional abundance predictions per spot
#
# NOTE: this MUST be run with Seurat v3 -- do not use v4!

# conda activate seurat3

# Input to provide from command line

# 1. Path to scRNA counts file (standard CytoSPACE input file format)
# 2. Path to cell type labels file (standard CytoSPACE input file format)
# 3. Path to ST data (standard CytoSPACE input file format)
# 4. Name of output file

# Run the script from command line with the specified inputs, i.e:
# Rscript get_cellfracs_seuratv3.R path_to_scrna_geneexpressionmatrix path_to_celltype_labels path_to_ST_geneexpressionmatrix name_of_output_file

##########################################################
library(data.table)
library(dplyr)
#require(devtools)
#install_version("Seurat",version="3.1.4",repos="http://cran.us.r-project.org")
library(Seurat)
library(ggplot2)
library(cowplot)

# Read inputs from command line
args = commandArgs(T)

fin_scrna = args[1]
fin_labels = args[2]
referencefile = args[3]
fn_out = args[4]

#####################################################
print("Reading in data")
expressions_noisy <- fread(fin_scrna, header = T, data.table = F)
rownames(expressions_noisy) = expressions_noisy[,1]
expressions_noisy = expressions_noisy[,-1]

annot <- read.delim(fin_labels)
annot = annot[match(colnames(expressions_noisy), annot$Cell.IDs),]
cell_index_vec = annot$CellType
names(cell_index_vec) = annot$Cell.IDs

reference <- fread(referencefile, header = T, data.table = F)
rownames(reference) = reference[,1]
reference = reference[,-1]
print("Done")
reference[is.na(reference)] <- 0
reference_col_sum <- as.matrix(colMeans(reference))
reference <- reference[,reference_col_sum > 0]
print("Creating Seurat object from ST data")
Seurat_reference <- CreateSeuratObject(counts = reference)
print("Done")

print("Performing SCTransform on ST data")
Seurat_reference <- SCTransform(Seurat_reference, assay = "RNA", verbose = FALSE, ncells=NULL) %>% RunPCA(verbose = FALSE)
print("Done")

expressions_noisy[is.na(expressions_noisy)] <- 0
expressions_row_sum <- as.matrix(rowMeans(expressions_noisy))
expressions_noisy <- expressions_noisy[expressions_row_sum > 0,] 

head(expressions_noisy)[,1:10]

print("Creating Seurat object from scRNA-seq data")
Seurat_expressions_noisy <- CreateSeuratObject(counts = expressions_noisy)
print("Done")
print("Performing SCTransform on scRNA-seq data")
Seurat_expressions_noisy <- SCTransform(Seurat_expressions_noisy, assay = "RNA", verbose = FALSE, ncells=NULL) %>% RunPCA()
print("Done")
print("Finding transfer anchors")
anchors <- FindTransferAnchors(reference = Seurat_expressions_noisy, query = Seurat_reference, normalization.method = "SCT")
print("Transfering data")
predictions.assay <- TransferData(anchorset = anchors, refdata = as.character(t(cell_index_vec)))
print("Done")

col_names <- sort(names(predictions.assay)[2:(dim(predictions.assay)[2]-1)])
col_sums <- vector()
for(i in 1:length(col_names)){
    col_sums[i] <- sum(select(predictions.assay,col_names[i]))
}
col_sum_total <- sum(col_sums)
col_sums_norm <- col_sums/col_sum_total

names(col_sums_norm) <- gsub('prediction.score.','',col_names)

names(col_sums_norm) <- gsub('.',' ',names(col_sums_norm),fixed=T)
df <- as.data.frame(col_sums_norm)
df$Index <- rownames(df)
colnames(df) <- c("Fraction","Index")
df <- df[,c("Index","Fraction")]

print("Writing estimated cell fractions to file")
write.table(t(df),file = fn_out, quote=F, sep="\t",row.names=T, col.names=F)
print("Done")

sessionInfo()

