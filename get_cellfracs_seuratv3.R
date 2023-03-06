suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("dplyr"))


#' @param fin_scrna Path to scRNA counts file (standard CytoSPACE input file format)
#' @param fin_labels Path to cell type labels file (standard CytoSPACE input file format)
#' @param fin_st Path to ST data (standard CytoSPACE input file format)
#' @param outdir Path to output directory
#' @param prefix Prefix of output file
#' @param st.norm.methods normalization methods for SRT data. CPM/SCT/log
#' @param sc.norm.methods normalization methods for scRNA-seq data. CPM/SCT/downsample/log
#' @param sc.max.umi Maximum number of UMIs, only used when performing downsampling normalization.
#' @import data.table dplyr Seurat
get_cellfracs_seuratv3 <- function(sc_path, ct_path, st_path, outdir){
    message(Sys.time(), " Load ST data")
    if(endsWith(st_path,'.csv')){
        st_delim <- ','
    } else{
        st_delim <- '\t'
    }
    st <- fread(st_path, sep = st_delim, header = TRUE, data.table = FALSE)
    rownames(st) = st[,1]; st = st[,-1]
    st[is.na(st)] <- 0
    st <- as.matrix(st)
    #st <- CreateSeuratObject(st) %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst") %>% ScaleData() %>% RunPCA(verbose = FALSE)
    st <- CreateSeuratObject(st) %>% SCTransform(verbose = FALSE, ncells = NULL) %>% RunPCA()
    message(Sys.time(), " Load scRNA data")
    if(endsWith(sc_path,'.csv')){
        sc_delim <- ','
    } else{
        sc_delim <- '\t'
    }
    scrna <- fread(sc_path, sep = sc_delim, header = TRUE, data.table = FALSE)
    rownames(scrna) = scrna[,1]; scrna = scrna[,-1]
    scrna[is.na(scrna)] <- 0
    num_cells <- dim(scrna)[2]
    scrna <- as.matrix(scrna)

    if(endsWith(ct_path,'.csv')){
        ct_delim <- ','
    } else{
        ct_delim <- '\t'
    }

    celltypes <- fread(ct_path, sep = ct_delim, header = TRUE, data.table = FALSE)
    rownames(celltypes) = celltypes[,1]
    colnames(celltypes)[2] <- "CellType"
    celltypes <- celltypes[colnames(scrna),]

    scrna <- CreateSeuratObject(counts = scrna, meta.data = celltypes) 
    Idents(scrna) <- scrna$CellType
    cell_types <- unique(scrna$CellType)
    num_cell_types <- length(cell_types)
    names(cell_types) <- make.names(cell_types)
    Idents(scrna) <- make.names(Idents(scrna))
    if(num_cells > 1000){
        scrna <- subset(scrna, downsample = 50)
    }
    #scrna <- NormalizeData(scrna) %>% FindVariableFeatures(selection.method = "vst") %>% ScaleData() %>% RunPCA(verbose = FALSE)
    scrna <- SCTransform(scrna, verbose = FALSE, ncells = NULL) %>% RunPCA()



    
    #print(unique(scrna$CellType))
    cell_index_vec <- Idents(scrna)
    message(Sys.time(), " Integration ")
    anchors <- FindTransferAnchors(reference = scrna, query = st, normalization.method = "SCT")
    predictions.assay <- TransferData(anchorset = anchors, refdata = cell_index_vec)
    colnames(predictions.assay) <- gsub('prediction.score.', '', colnames(predictions.assay))
    predictions.assay <- predictions.assay[, 2:(ncol(predictions.assay)-1)]
    col_sums <- colSums(predictions.assay)
    cellfrac <- col_sums / sum(col_sums)
    cellfrac <- data.frame(Index = names(cellfrac), Fraction = cellfrac)
    cellfrac$Index <- cell_types[cellfrac$Index]

    write.table(predictions.assay, paste0(outdir, "/Seurat_cellfracs.txt"), quote = FALSE, sep = "\t")
    write.table(t(cellfrac), paste0(outdir, "/Seurat_weights.txt"), quote = FALSE, sep = "\t", col.names = FALSE)
} 


parser <- ArgumentParser()
# specify options 
parser$add_argument("--scrna-path", type="character", default="",
                    help="Path to single cell RNA-seq data, with rows as genes and columns as single cells")
parser$add_argument("--ct-path", type="character", default="",
                    help="Path to cell type annotation file, with first column showing single cell names, and second column showing cell type annotations")
parser$add_argument("--st-path", type="character", default="",
                    help="Path to Spatially-Resolved Transcriptomics (SRT) data, with rows as genes and columns as spots")
parser$add_argument("--outdir", type="character", default="./", help="Output directory")
parser$add_argument("--prefix", type="character", default="CytoSPACE_input.", help="Prefix for output files")


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()
get_cellfracs_seuratv3(sc_path = args$scrna_path, ct_path = args$ct_path, 
                       st_path = args$st_path, outdir = args$outdir)
