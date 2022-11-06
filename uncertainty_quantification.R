#' Evaluating the confidence of CytoSPACE alignment of single cells to spatial spots.
#' @examples 
#' Run the script from command line with the specified inputs, i.e:
#' Rscript uncertainty_quantification.R path/to/ST_geneexpressionmatrix path/to/scrna_geneexpressionmatrix path/to/assigned_locations.csv name_of_output_file(optional)
#'
#' @param 
#' Input to provide from command line
#' 1. Path to ST data (standard CytoSPACE input file format)
#' 2. Path to scRNA counts file (standard CytoSPACE input file format)
#' 3. Path to the assigned_locations.csv file (standard CytoSPACE output)
#' 4. Name of output file
#' 
#' @details
#' Pseudo-code
#' 1. Calculate the top 50 marker genes by –log10 adjusted p-value per annotated cell type in the scRNA-seq query dataset.
#' 2. For each cell type i mapped to the ST sample:
#'     a. Use the scRNA-seq query dataset to reconstitute pseudo-bulk spot transcriptomes that mirror single-cell assignments. Repeat for 50 randomly selected spots containing at least one cell of cell type i and 50 randomly selected spots lacking at least one cell of cell type i.
#'     b. Train a support vector machine (SVM) model to distinguish the two groups from the previous stem using the top 50 marker genes of cell type i
#'     c. Apply the SVM model to calculate the probability of class assignment, termed a confidence score, for each spot in the original ST sample.
#'     d. For each mapped cell of type i, retrieve its spot-specific confidence score.
#'

library(Seurat) #tested with v4.0.1
library(data.table) #tested with v1.14.0
library(e1071) #tested with v1.7.8

# Read inputs from command line
args = commandArgs(T)

#========================Input files============================
#ST dataset (count matrix)
stfile <- args[1]
# stfile <- "CytoSPACE_example_breast_cancer/brca_STdata_GEP.txt"
#scRNA-seq dataset (count matrix)
scfile <- args[2]
# scfile <- "CytoSPACE_example_breast_cancer/brca_scRNA_GEP.txt"
#CytoSPACE output
assignfile <- args[3]
# assignfile <- "assigned_locations.csv"
# output file name
outputname <- gsub(".csv", "wConfidenceScores.csv", assignfile)
if(length(args)>3) outputname <- args[4]

#=========================Parameters===========================

groupsize <- 50 #number of pseudobulks per group
nummarkers <- 50 #number of marker genes to use per cell type (maximum)
seednum <- 1 #default seed


#=================Read in files and preprocess==================

#ST
print(">Reading in ST dataset...")
st <- fread(stfile, header=T)
st <- data.frame(st)
rownames(st) <- make.names(st[,1],unique=T)
st <- st[,-1]

#SC
print(">Reading in scRNA-seq dataset...")
sc <- fread(scfile, header=T)
sc <- data.frame(sc)
rownames(sc) <- make.names(sc[,1],unique=T)
sc <- sc[,-1]

inter <- intersect(rownames(st),rownames(sc))
st <- st[inter,]
sc <- sc[inter,]

#CytoSPACE output
assign <- read.csv(assignfile)

cells <- unique(assign$CellType) #unique mapped cell types

sc <- sc[,assign$OriginalCID]
colnames(sc) <- assign$UniqueCID

spots <- unique(assign$SpotID)
st <- st[,which(colnames(st) %in% spots)]

sclabels <- assign$CellType

#Normalize data and identify DEGs

print(">Processing data with Seurat...")
x <- CreateSeuratObject(sc)
Idents(x) <- sclabels
x <- NormalizeData(x)
x <- ScaleData(x)
sc  <- x[["RNA"]]@scale.data
sc <- data.matrix(sc)
set.seed(seednum)
if(ncol(sc)>500) x <- subset(x, downsample=50) #downsample if large before finding markers
markers <- FindAllMarkers(x)
markers <- markers[which(markers$avg_log2FC>0),]

x <- CreateSeuratObject(st)
x <- NormalizeData(x)
x <- ScaleData(x)
st  <- x[["RNA"]]@scale.data
st <- data.matrix(st)

#=================Model training and application==================

conf <- matrix(rep(0,length(cells)*ncol(st)),nrow=ncol(st), byrow=T) #save confidence scores for each cell type

#iterate through each mapped cell type in the scRNA-seq dataset
print(">SVM training and application...")

for(index in 1:length(cells)){

    print(cells[index])
    
    gt <- rep(0,ncol(st))
    gt[which(colnames(st) %in% assign$SpotID[which(assign$CellType==cells[index])])] <- 1
    names(gt) <- colnames(st)

    #marker genes
    mk <- markers$gene[which(markers$cluster==cells[index])]
    mk <-  mk[head(which(mk %in% rownames(st)), nummarkers)] #select top markers
    
    #if at least 5 marker genes are available
    
    if(length(mk) >= 5){

        sc2 <- sc

        #Create pseudobulks (scr1) for spots containing at least 1 cell of the current cell type
        
        spots <- unique(names(gt)[which(gt==1)])
        set.seed(seednum)
        if(length(spots)>groupsize) spots <- s1 <- sample(spots, groupsize)
        if(length(spots)<groupsize) spots <- s1 <- sample(spots, groupsize, replace = T)
        
        for(i in 1:length(spots)){
            if(length(which(assign$SpotID == spots[i]))>1){
                d <- apply(sc2[,which(colnames(sc) %in% assign$UniqueCID[which(assign$SpotID == spots[i])])], 1, mean)
            }else{
                d <- sc2[,which(colnames(sc) %in% assign$UniqueCID[which(assign$SpotID == spots[i])])]
            }
            if(i == 1) scr <- d
            else{scr <- cbind(scr,d)}
        }

        scr1 <- scr
        colnames(scr1) <- spots
        
        #Create pseudobulks (scr2) for spots not containing at least 1 cell of the current cell type
        
        spots <- unique(names(gt)[which(gt==0)])
        set.seed(seednum)
        if(length(spots)>groupsize)  spots <- s2 <- sample(spots, groupsize)
        if(length(spots)<groupsize) spots <- s2 <- sample(spots, groupsize, replace = T)

        for(i in 1:length(spots)){
            if(length(which(assign$SpotID == spots[i]))>1){
                d <- apply(sc2[,which(colnames(sc) %in% assign$UniqueCID[which(assign$SpotID == spots[i])])], 1, mean)
            }else{
                d <- sc2[,which(colnames(sc) %in% assign$UniqueCID[which(assign$SpotID == spots[i])])]
            }
            if(i == 1) scr <- d
            else{scr <- cbind(scr,d)}
        }

        scr2 <- scr
        colnames(scr2) <- spots

        #Model training/validation

        #Train SVM model on pseudobulks
        train <- t(cbind(scr1,scr2)[mk,])
        df <- data.frame(class=as.factor(c(rep("Positive",groupsize),rep("Negative",groupsize))), train)
        model <- svm(class ~ ., data = df, probability=TRUE)

        #Application of SVM model to ST dataset
        test <- data.frame(t(st[mk,]))
        colnames(test) <- colnames(t(st[mk,]))
        
        pred <- predict(model, test, probability=TRUE)
        out <- attr(pred, "probabilities")[,1]
        
        conf[,index] <- out
    
    }

}

#Compile output and write to disk

colnames(conf) <- cells
rownames(conf) <- colnames(st)

assign_all <- assign
assign_all <- assign_all[which(assign_all$SpotID %in% colnames(st)),]
confid <- sapply(1:nrow(assign_all), function(i) conf[assign_all$SpotID[i],assign_all$CellType[i]])
assign_all <- cbind(assign_all, confid)
colnames(assign_all)[ncol(assign_all)] <- "Confidence" #append confidence scores to CytoSPACE output

write.csv(assign_all, file = outputname, row.names=F, quote=F) #write to disk

print(">Done")
