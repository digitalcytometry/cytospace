<img src="https://github.com/digitalcytometry/cytospace/blob/main/images/CytoSPACE_logo.jpeg" width="300"> 

# CytoSPACE: Robust alignment of single-cell and spatial transcriptomes

**CytoSPACE** is a novel computational strategy for assigning single-cell transcriptomes to in situ spatial transcriptomics (ST) data. Our method solves single cell/spot assignment by minimizing a correlation-based cost function through a linear programming-based optimization routine. 

<img src="https://github.com/digitalcytometry/cytospace/blob/main/images/CytoSPACE_overview.png" width="900"> 

The key innovations of our method are:

- In contrast to conventional methods, CytoSPACE dissects spatial organizations of cells in a given tissue at single cell level.

- Since our method maps single cells from scRNA-sequencing data, in which larger numbers of genes are sequenced per each cell compared to available spatial transcriptomics technology, our method improves the gene coverage of a recontructed tissue significantly.

## Installation instructions
1. Install <a href="https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html" target="_blank">Miniconda</a> if not already available.

2. Clone this repository (`git clone`)

3. Navigate to `cytospace` directory:
```bash
cd cytospace
```
4. Create a conda environment with the required dependencies:
```bash
conda env create -f environment_withoutlapjv.yml
```
5. Activate the `cytospace` environment you just created:
```bash
conda activate cytospace
``` 

6. Install CytoSPACE by executing:
```bash
pip install .
``` 

7. (Optional) Install package `lapjv` by executing:
```bash
pip install lapjv
```
We highly recommend you install this package, which provides a fast implementation of the core optimization algorithm within CytoSPACE. However, some systems may not accommodate it as it requires CPU support for AVX2 instructions. To determine if your system supports this package, it is generally easiest to simply attempt to install it as above. If you run into an error, it is likely your system does not support it, and you can simply use one of the other options we have provided. See __Solver options__ below for details. 

## File format
CytoSPACE requires 5 files as input. All files should be provided in tab or comma-delimited tabular input format (saved as .txt or .csv, respectively) with no double quotations. Further formatting details for each input file are specified below:

1. __A scRNA-seq gene expression file:__
- The matrix must be genes (rows) by cells (columns).
- The first row must contain the single cell IDs and the first column must contain the gene names.
- The gene expression data should be represented as non-normalized counts. 
<img src="https://github.com/digitalcytometry/cytospace/blob/main/images/scRNAfile.png" width="800"> 

2. __A cell type label file:__
- Cell type labels corresponding to the single cell IDs in the scRNA-seq gene expression matrix. 
- The table should contain two columns, where column 1 contains the single cell IDs corresponding to the columns of the scRNA-seq matrix and column 2 contains the corresponding cell type labels.
- The columns must have a header. 
<img src="https://github.com/digitalcytometry/cytospace/blob/main/images/celllabelfile.png" width="250"> 

3. __A spatial transcriptomics (ST) gene expression file:__
- The matrix must be genes (rows) by ST spots (columns).
- The first row must contain the ST spot IDs and the first column must contain the gene names.
- The gene expression data should be represented as non-normalized counts. 
<img src="https://github.com/digitalcytometry/cytospace/blob/main/images/STdatafile.png" width="800"> 

4. __A spatial transcriptomics coordinates file:__
- A table consisting of 3 columns, where the first column contains the ST spot IDs corresponding to the columns of the ST gene expression matrix, and column 2 and 3 contain the X and Y coordinates of the spatial transcriptomics data, respectively. 
- The columns must have a header. 
<img src="https://github.com/digitalcytometry/cytospace/blob/main/images/STcoordfile.png" width="250"> 

5. __A file with cell type fraction estimates, obtained from the `R` script `get_cellfracs_seuratv3.R`.__ 
- A table consisting of 2 rows, where the first row is the cell type labels, and the second row is the cell fractions of each cell type represented as proportions between 0 and 1. The first column is the row names. 
- For further details on running `get_cellfracs_seuratv3.R`, see section "__Preprocessing__" below.

## Preprocessing
To account for the disparity between scRNA-seq and ST data in the number of cells per cell type, the fractional composition of each cell type in the ST tissue needs to be provided as input to CytoSPACE. This is determined using an external deconvolution tool, such as Spatial Seurat, CIBERSORTx, or SPOTlight. We have included Spatial Seurat in our benchmarking, and provide here a script to obtain the cell type fractions using this approach.

Run the script `get_cellfracs_seuratv3.R` from command line with the following inputs:
1. Path to scRNA counts file (same scRNA-seq gene expression matrix input file format as specified in __File format__ section point 1)
2. Path to cell type labels file (same cell type label input file format as specified above in __File format__ section point  2)
3. Path to ST data (same ST gene expression matrix input file format as specified above in __File format__ section point  3)
4. Name of output file

For example:
```bash
Rscript /path/to/get_cellfracs_seuratv3.R melanoma_scRNA_GEP.txt melanoma_scRNA_celllabels.txt melanoma_STdata_slide1_GEP.txt melanoma_cell_fraction_estimates.txt
```
### Important, please note:
1. While `cytospace` can be run from any path and folder, the path to `get_cellfracs_seuratv3.R` must be specified in the command. 
2. We use `Seurat v3` for estimating cell fractions, and is installed as part of the CytoSPACE environment. If you want to run other analyses using more recent versions of Seurat after running CytoSPACE, for example Seurvat v4, make sure to first __deactivate the CytoSPACE environment__ once you are done running CytoSPACE. This is done using the command `deactivate cytospace`.

## Running CytoSPACE
CytoSPACE can be called from the command line from any folder using `cytospace`. Examples on how to run CytoSPACE are provided in the section "Example datasets for running CytoSPACE" below.

A typical CytoSPACE run with default settings would look like this: 
 ```bash
 cytospace --scRNA-path /path/to/scRNA_geneexpression
    --cell-type-path /path/to/scRNA_celllabels
    --st-path /path/to/ST_geneexpression
    --coordinates-path /path/to/ST_coordinates
    --cell-type-fraction-estimation-path path/to/cellfracestimates
```
Or with more condensed parameter names: 
 ```bash
 cytospace -sp /path/to/scRNA_geneexpression
    -ctp /path/to/scRNA_celllabels
    -stp /path/to/ST_geneexpression
    -cp /path/to/ST_coordinates
    -ctfep path/to/cellfracestimates
```
To see a list of variables and default values for running CytoSPACE, you can call `cytospace` from the command line along with the `-h` or 
`--help` flag, i.e., `cytospace -h`.

### Other ways CytoSPACE can be run:
1. You can call the `cytospace.py` script directly with python:
 `python cytospace/cytospace.py`
 
2. You can import methods or functions from `CytoSPACE` in python and modify/create your own 
    pipeline. For example:
```python
from cytospace import cytospace

for mean_cell_numbers in [5, 10, 20]:
    cytospace(..., mean_cell_numbers=mean_cell_numbers)
```

## Example datasets for running CytoSPACE
For users to test CytoSPACE, we have included files for two example runs:
1. A melanoma scRNA-seq atlas by Tirosh et al (<a href="https://www.science.org/doi/10.1126/science.aad0501?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed" target="_blank">Science, 2016</a>), and a melanoma specimen profiled by the legacy ST platform (Thrane et al, <a href="https://aacrjournals.org/cancerres/article/78/20/5970/631815/Spatially-Resolved-Transcriptomics-Enables" target="_blank">Cancer Research, 2018</a>). This example is very quick to run, and a good test case to make sure that CytoSPACE is running as expected.
2. A HER2+ breast cancer scRNA-seq atlas by Wu et al (<a href="https://www.nature.com/articles/s41588-021-00911-1" target="_blank">Nature Genetics, 2021</a>) and a HER2+ breast cancer FFPE specimen profiled by the Visium platform (<a href="https://www.10xgenomics.com/resources/datasets/human-breast-cancer-ductal-carcinoma-in-situ-invasive-carcinoma-ffpe-1-standard-1-3-0" target="_blank">10x Genomics</a>)

### Download example datasets
Zip files containing the example datasets can be downloaded from the following links:
1. <a href="https://drive.google.com/file/d/1oHe4UP2K0kQS9gNFVvZtjqtJyJK_VeBp/view?usp=sharing" target="_blank">Melanoma</a> 
2. <a href="https://drive.google.com/file/d/1vAqszYk3-B2vgwkSFMprsUcRBFr-lS2f/view?usp=sharing" target="_blank">Breast cancer</a>

To download from the command line using `gdown`:
1. Melanoma
   ```bash
   gdown --fuzzy https://drive.google.com/file/d/1oHe4UP2K0kQS9gNFVvZtjqtJyJK_VeBp/view?usp=sharing
   unzip CytoSPACE_example_melanoma.zip
   ```
2. Breast cancer
   ```bash
   gdown --fuzzy https://drive.google.com/file/d/1vAqszYk3-B2vgwkSFMprsUcRBFr-lS2f/view?usp=sharing
   unzip CytoSPACE_example_breast_cancer.zip
   ```
   
   
### Commands for running example analyses:
Once the example files are downloaded, the commands below can be run from the folders where the example datasets are located:
 ```bash
cytospace -sp melanoma_scRNA_GEP.txt -ctp melanoma_scRNA_celllabels.txt -stp melanoma_STdata_slide1_GEP.txt -cp melanoma_STdata_slide1_coordinates.txt -ctfep melanoma_cell_fraction_estimates.txt
```

```bash
cytospace -sp brca_scRNA_GEP.txt -ctp brca_scRNA_celllabels.txt -stp brca_STdata_slide1_GEP.txt -cp brca_STdata_slide1_coordinates.txt -ctfep brca_cell_fraction_estimates.txt
```
## CytoSPACE output files
The main output from a CytoSPACE run is the file named `assigned_locations.csv`, which provides the ST spots to which the single cells have been assigned.

The CytoSPACE results are visualized in heatmaps saved as `plot_cell_type_locations.pdf` showing the distribution of single cells across ST spots for each cell type. 

## Authors
CytoSPACE was developed by

* Milad R. Vahid (miladrv)
* Erin L. Brown (erinlbrown)
* Chloé B. Steen (cbsteen)
* Aaron M. Newman (aaronmnewman)

## Licence
CytoSPACE is licensed under the GNU GPL, version 3 or (at your option) any
later version.
CytoSPACE is Copyright (2022-) by the authors.

## Citation
If you are using CytoSPACE for scientific work, please cite it by:

    TODO
