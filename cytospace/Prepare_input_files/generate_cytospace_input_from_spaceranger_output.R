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
##########################################################
library(hdf5r)
library(rhdf5)
library(Seurat)
library(Matrix)

# Read inputs from command line
args = commandArgs(T)

# fin_ST = args[1]
# fn_out = args[2]

fin_ST <- '/Users/i0461476/Desktop/OpenTargets/Spatial_data/data/'
fn_out <- '/Users/i0461476/Desktop/OpenTargets/Spatial_data/data/'

ST_data <- Load10X_Spatial(fin_ST)
ST_expressions <- as.matrix(ST_data@assays$Spatial@counts)
ST_expressions <- cbind(rownames(ST_expressions), ST_expressions)
colnames(ST_expressions)[1] <- 'GENES'
coordinates_raw <- ST_data@images$slice1@coordinates
coordinates <- coordinates_raw[,-1]
coordinates <- cbind(rownames(coordinates), coordinates)
colnames(coordinates)[1] <- 'Spot ID'

print("Writing output to file")
write.csv(coordinates, file = paste(fn_out, 'Coordinates.csv'), row.names = F)
write.csv(ST_expressions, file = paste(fn_out, 'ST_data.csv'), row.names = F)
print("Done")

sessionInfo()