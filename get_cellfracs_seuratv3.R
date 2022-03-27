
##########################################################
# README
##########################################################
# Script to compute overall cell type fractions for ST
# sample given scRNA set with cell type annotations as
# well as fractional abundance predictions per spot
#
# NOTE: We recommend that this be run with Seurat v3
#
# Input to provide from command line
#
# 1. Path to scRNA counts file (standard CytoSPACE input file format)
# 2. Path to cell type labels file (standard CytoSPACE input file format)
# 3. Path to ST data (standard CytoSPACE input file format)
# 4. Name of output file
#
# Run the script from command line with the specified inputs, i.e:
# Rscript get_cellfracs_seuratv3.R path_to_scrna_geneexpressionmatrix path_to_celltype_labels path_to_ST_geneexpressionmatrix name_of_output_file

##########################################################
library(data.table)
library(dplyr)
library(Seurat)

# Read inputs from command line
args = commandArgs(T)

fin_scrna = args[1]
fin_labels = args[2]
fin_st = args[3]
fn_out = args[4]

#####################################################
print("Reading in data")
scrna <- fread(fin_scrna, header = T, data.table = F)
rownames(scrna) = scrna[,1]
scrna = scrna[,-1]

annot <- read.delim(fin_labels)
cell_ids <- annot[,1]
cell_types <- annot[,2]
annot = annot[match(colnames(scrna), cell_ids),]
cell_index_vec = cell_types
names(cell_index_vec) = cell_ids

st <- fread(fin_st, header = T, data.table = F)
rownames(st) = st[,1]
st = st[,-1]
print("Done")
st[is.na(st)] <- 0
st_col_sum <- as.matrix(colMeans(st))
st <- st[,st_col_sum > 0]
print("Creating Seurat object from ST data")
st <- CreateSeuratObject(counts = st)
print("Done")

print("Performing SCTransform on ST data")
st <- SCTransform(st, assay = "RNA", verbose = FALSE, ncells=NULL) %>% RunPCA(verbose = FALSE)
print("Done")

scrna[is.na(scrna)] <- 0
scrna_row_sum <- as.matrix(rowMeans(scrna))
scrna <- scrna[scrna_row_sum > 0,] 

head(scrna)[,1:10]

print("Creating Seurat object from scRNA-seq data")
scrna <- CreateSeuratObject(counts = scrna)
print("Done")
print("Performing SCTransform on scRNA-seq data")
scrna <- SCTransform(scrna, assay = "RNA", verbose = FALSE, ncells=NULL) %>% RunPCA()
print("Done")
print("Finding transfer anchors")
anchors <- FindTransferAnchors(reference = scrna, query = st, normalization.method = "SCT")
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

