<p align="center">
  <img width="300" src="images/CytoSPACE_logo.jpeg">
</p>

<h1> <p align="center">
    Robust and rapid alignment of single-cell and spatial transcriptomes
</p> </h1>

**CytoSPACE** is a novel computational tool for assigning single-cell transcriptomes to in situ spatial transcriptomics (ST) data. Our method solves single cell/spot assignment by minimizing a correlation-based cost function through a shortest augmenting path optimization routine. 

<p align="center">
  <img src="images/CytoSPACE_overview.png" width="900"> 
</p>

The key innovations of our method are:

- Unlike conventional methods which calculate cell type decompositions by spot, CytoSPACE yields a reconstructed tissue specimen with both high gene coverage and spatially-resolved scRNA-seq data suitable for downstream analysis.
- CytoSPACE is highly robust to noise and returns globally optimal cell-to-spot assignments. (See the [paper](https://www.nature.com/articles/s41587-023-01697-9) for full details.)
- Unlike other methods which generally operate on pre-selected marker genes or on a shared embedding space (the latter of which can erase true biological variation), CytoSPACE uses the full transcriptome without the need for batch correction, helping it retain sensitivity to subtle cell states.
- CytoSPACE is quick and simple to execute. It runs in minutes even with a single CPU on a personal laptop and requires no hyperparameter tuning or gene/feature selection.

CytoSPACE is available through a web interface at <a href="https://cytospace.stanford.edu/">cytospace.stanford.edu</a>, which enables users to run CytoSPACE with default settings without downloading the source code.

## Installation instructions (5-10 minutes)

<details><summary>Expand section</summary>

1. Install <a href="https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html" target="_blank">Miniconda</a> if not already available.

2. Clone this repository:
```bash
  git clone https://github.com/digitalcytometry/cytospace
```

3. Navigate to `cytospace` directory:
```bash
  cd cytospace
```

4. (5-10 minutes) Create a conda environment with the required dependencies:
```bash
  conda env create -f environment.yml
```

5. Activate the `cytospace` environment you just created:
```bash
  conda activate cytospace
``` 

6. (~30 seconds) Install CytoSPACE by executing:
```bash
  pip install .
``` 

7. (Recommended, ~1 minute) Install package `lapjv` by executing:
```bash
  pip install lapjv==1.3.14
```
We highly recommend you install this package, which provides a fast implementation of the default core optimization algorithm within CytoSPACE.

__Step 7 may result in an error depending on your system__, as the package requires CPU support for AVX2 instructions that not all systems support. To determine if your system supports this package, it is generally easiest to simply attempt to install it as above. If it installs without problems, your system will support it! If you run into an error, it is likely that your system does not support it, and you can simply use one of the other options we have provided by specifying the `--solver-method` (`-sm`) parameter in your call to CytoSPACE. See [__Running CytoSPACE__](#running-cytospace) - __Choosing a solver__ for details.

Please note that, if the package successfully installs but you receive an "illegal instruction" error while running CytoSPACE, you may be able to build the package instead with the following command:
```bash
   pip3 install git+https://github.com/src-d/lapjv
```
For more information, see the <a href="https://pypi.org/project/lapjv/" target="_blank">lapjv documentation page</a>. 
</details>

## Input files

<details><summary>Expand section</summary>

By default, CytoSPACE requires 4 files as input. All files should be provided in a tab-delimited tabular input format (.txt) with no double quotations. Further formatting details for each input file are specified below. At the end of this section, we also provide instructions on using scripts to generate input files from Seurat objects.

1. __A scRNA-seq gene expression file:__
- The matrix must be genes (rows) by cells (columns).
- The first row must contain the single cell IDs and the first column must contain the gene names.
- The first column (gene names) must have a header.
- The gene expression data should be represented as non-normalized counts. 
- All instances of duplicate gene names will be dropped at runtime.
<p align="center">
  <img src="images/scRNAfile.png" width="800"> 
</p>

2. __A cell type label file:__
- Cell type labels corresponding to the single cell IDs in the scRNA-seq gene expression matrix. 
- Cell type label strings should not have special characters. 
- The table should contain two columns, where column 1 contains the single cell IDs corresponding to the columns of the scRNA-seq matrix and column 2 contains the corresponding cell type labels.
- The columns must have a header. 
<p align="center">
  <img src="images/celllabelfile.png" width="300"> 
</p>

3. __A spatial transcriptomics (ST) gene expression file:__
- The matrix must be genes (rows) by ST spots (columns).
- The first row must contain the ST spot IDs and the first column must contain the gene names.
- The first column (gene names) must have a header.
- The gene expression data should be represented as non-normalized counts. 
- All instances of duplicate gene names will be dropped at runtime.
<p align="center">
  <img src="images/STdatafile.png" width="800"> 
</p>

4. __A spatial transcriptomics coordinates file:__
- A table consisting of 3 columns, where the first column contains the ST spot IDs corresponding to the columns of the ST gene expression matrix and columns 2 and 3 contain the row and column indices of each ST spot, respectively. 
- The columns must have a header. 
<p align="center">
  <img src="images/STcoordfile.png" width="300"> 
</p>

### From Space Ranger outputs
If the users are starting from Space Ranger outputs, they can provide the ST input files as a single tar.gz, __in place of__ both (3) gene expression and (4) coordinates. If a Space Ranger output is specified, CytoSPACE will automatically attempt to unzip the provided tarball and load the correponding ST expression and coordinates data.<br>
If there are duplicate gene names in the gene expression file, they will be made unique as is done in Seurat.

The tarball should only include the following:
- A single H5 file (extension .h5) containing the ST gene expression
- A single subdirectory containing the image data

With the above items inside a directory named `spaceranger_input`, a tarball can be generated using the following command:
```bash
  tar -cvzf sr_input.tar.gz spaceranger_input
```
Or more generally:
```bash
  tar -cvzf [name_of_tarball] [name_of_directory]
```

An example file tree for an unzipped tarball is shown below on the left. If downloading from the public 10X Visium data, users can download the files shown below on the right.
<p align="center">
  <img src="images/VisiumTar.png" width="300"> <img src="images/Visium.png" width="300">
</p>

Starting with CytoSPACE v1.1.0, the tool supports Space Ranger v2.0.0+ outputs. This is done by detecting files that have updated naming conventions (namely, the `tissue_positions.csv` file), then saving a copy of those files to the same directory with a naming convention that is compatible with scipy's functionality to read Visium datasets.

### Using sparse matrices for gene expression
Starting with CytoSPACE v1.0.5, users may also provide the scRNA or ST gene expression as sparse matrices instead of tab- or comma-delimited files.<br> 
This will require a very specific set of file names in hopes of avoiding issues with parsing, and we recommend that users use the R helper scripts under `Prepare_input_files` to generate these input.  Please see the subsections below for further information about these helper scripts.

If providing input as sparse matrices, you will need three files under the same directory to represent one expression matrix: `[expression].mtx`, `genes.tsv`, and `cells.tsv`. The `.mtx` file will list the numeric values in sparse matrix format, and can also be in a compressed `[expression].mtx.gz` format. The filename inside the brackets (`expression`) may vary. The `genes.tsv` and `cells.tsv` will be the corresponding gene names and cell (or spot/barcode) names for the matrix, located in the same directory as the matrix file. 

Starting CytoSPACE v1.1.0, the tool also supports alternative naming conventions (`features.tsv` and `barcodes.tsv`), that are according to common Space Ranger and Cell Ranger sparse matrix outputs. These files support `.csv` extension, and compressed `.gz` formats as well. The files should have one entry on each line (no headers). If there are multiple columns in any of the files, CytoSPACE will take the first column as the gene/cell names. Please note that for common Cell Ranger and Space Ranger sparse matrix inputs, the first column of the gene file usually has the Ensembl IDs, while the gene symbols are in the second column, so you would need to change the order of the columns if your other input modality contains gene symbols for the gene names.

The `[expression].mtx` file should be supplied as the argument to `--st-path` or `--scRNA-path`, in which case CytoSPACE will automatically try to locate the corresponding gene and cell files from the same directory.

Please ensure the gene file includes gene names with the same naming convention (e.g., Ensembl gene ID, HGNC gene symbol, etc.) as those in the provided ST data.

<details><summary><b>Preparing input files from Seurat objects</b></summary>

If you have data in the form of Seurat objects, you can generate files formatted for CytoSPACE input via helper functions we have provided in the `R` script `generate_cytospace_from_seurat_object.R` in `cytospace/Prepare_input_files`. To use these helper functions, first import them from `generate_cytospace_from_seurat_object.R` by including 
```bash
  source('/path/to/generate_cytospace_from_seurat_object.R')
```
in your R script. 

### From scRNA-seq Seurat object
For producing CytoSPACE inputs from scRNA Seurat objects, we provide the function `generate_cytospace_from_scRNA_seurat_object` which may be called as
```bash
  generate_cytospace_from_scRNA_seurat_object(scRNA_Seurat_Object, dir_out='.', fout_prefix='', write_sparse=FALSE, rna_assay='RNA')
```
within your R script.<br>
`scRNA_Seurat_Object` (required) : input Seurat object<br>
`dir_out` (optional, default is working directory) : the path to the output directory to store the results<br>
`fout_prefix` (optional, default is none) : a prefix to add to output file names, which otherwise are generated as `scRNA_data.txt` and `cell_type_labels.txt`<br>
`write_sparse` (optional, default is FALSE) : whether to save the expression data in sparse matrix format<br>
`rna_assay` (optional, default is `RNA`) : which assay to take the count matrix from<br>
Please note that `Idents(scRNA_Seurat_Object)` must be set to include cell types.


### From Spatial Seurat object
For producing CytoSPACE inputs from ST Seurat objects, we provide the function `generate_cytospace_from_ST_seurat_object` which may be called as
```bash
  generate_cytospace_from_ST_seurat_object(ST_Seurat_Object, dir_out='.', fout_prefix='', write_sparse=FALSE, slice='slice1')
```
within your R script.<br>
`ST_Seurat_Object` (required) : input Seurat object<br>
`dir_out` (optional, default is working directory) : the path to the output directory to store the results<br>
`fout_prefix` (optional, default is none) : a prefix to add to output file names, which otherwise are generated as `ST_data.txt` and `Coordinates.txt`<br>
`write_sparse` (optional, default is FALSE) : whether to save the expression data in sparse matrix format<br>
`slice` (optional, default is `slice1`) provides the name of your slice as stored in your Seurat object.
</details>
</details>

## Running CytoSPACE

<details><summary>Expand section</summary>

After activating the `cytospace` conda environment via `conda activate cytospace`, CytoSPACE can be called from the command line from any folder using `cytospace`. More examples on how to run CytoSPACE are provided in the [__Example dataset for running CytoSPACE__](#example-dataset-for-running-cytospace) section below.

A typical CytoSPACE run with default settings would look like this: 
 ```bash
 cytospace \
    --scRNA-path /path/to/scRNA_geneexpression \
    --cell-type-path /path/to/scRNA_celllabels \
    --st-path /path/to/ST_geneexpression \
    --coordinates-path /path/to/ST_coordinates
```
Or with more condensed parameter names: 
 ```bash
 cytospace \
    -sp /path/to/scRNA_geneexpression \
    -ctp /path/to/scRNA_celllabels \
    -stp /path/to/ST_geneexpression \
    -cp /path/to/ST_coordinates
```

Alternatively, if starting from a Space Ranger output, the command may look like this:
```bash
 cytospace \
    --scRNA-path /path/to/scRNA_geneexpression \
    --cell-type-path /path/to/scRNA_celllabels \
    --spaceranger-path /path/to/spaceranger_output.tar.gz
```
```bash
 cytospace -sp /path/to/scRNA_geneexpression \
    -ctp /path/to/scRNA_celllabels \
    -srp /path/to/spaceranger_output.tar.gz
```

For full usage details with additional options, see the [__Extended usage details__](#extended-usage-details) section below. 

### Choosing a solver
CytoSPACE provides three solver options. In short, we recommend using the default option `lapjv` if your system supports AVX2 (i.e., if you were able to successfully install it with `pip install lapjv==1.3.14` on Step 7 of the installation) and `lap_CSPR` otherwise. No options are required to use the default solver `lapjv`. To use `lap_CSPR` instead, pass the argument `-sm lap_CSPR` to your `cytospace` call. For full solver details, see the [__CytoSPACE solver options__](#cytospace-solver-options) section below.

<details><summary><b>Other ways CytoSPACE can be run</b></summary>
 
- You can import methods or functions from `CytoSPACE` in python and modify/create your own 
    pipeline. For example:
```python
  from cytospace import cytospace

  for mean_cell_numbers in [5, 10, 20]:
      cytospace.main_cytospace(..., mean_cell_numbers=mean_cell_numbers)
```
</details>
</details>

## CytoSPACE outputs

<details><summary>Expand section</summary>

CytoSPACE will produce six output files by default.
1. ```cell_type_assignments_by_spot.pdf```<br>
Heatmaps of cell type assignments within the ST sample. Along with a plot showing the total number of cells mapped to each spot, these show the spatial distribution of cell type assignments. Color bars indicate the number of cells of the respective cell type inferred per spot.
2. ```cell_type_assignments_by_spot_jitter.pdf```<br>
A single scatterplot showing all assigned cells by their spot location. Each cell is colored based on its cell type.
3. ```assigned_locations.csv```<br>
This file will provide the assigned locations of each single cell mapped to ST spots. As some cells may be mapped to multiple locations depending on the size of the input scRNA-seq set, new cell IDs (`UniqueCID`) are assigned to each cell and given in the first column. The second column includes original cell IDs (`OriginalCID`); the third column includes corresponding cell types (`CellType`); the fourth column includes assigned spot IDs (`SpotID`); and the fifth and sixth columns respectively include  `row` and `column` indices, or xy-coordinates such as `X` and `Y` if provided in the initial coordinates file, of the corresponding spots.
4. ```assigned_expression```, a directory with `barcodes.tsv`, `genes.tsv`, and `matrix.mtx` (CytoSPACE v1.0.6+)<br>
This represents the gene expression of the resulting assignments, with rows as genes and columns as cells (`UniqueCID` of `assigned_locations.csv`), as recovered from the original input scRNA matrix. The expression data can be read in for downstream analyses using functions such as Seurat `Read10X()` or SciPy `io.mmread()`. As with `assigned_locations.csv`, please note that there may be cells that are assigned to multiple locations and therefore appear multiple times in the expression matrix (under different `UniqueCID`s).<br>
For compatibility with the default parameters of Seurat `Read10X()`, `genes.tsv` lists the gene names in its second column, with the first column filled with `NA`s.
5. ```cell_type_assignments_by_spot.csv```<br>
This file gives the raw number of cells of each cell type per spot by `SpotID` as well as the total number of cells assigned to that spot.
6. ```fractional_abundances_by_spot.csv```<br>
This file gives the fractional abundance of cell types assigned to each spot by `SpotID`.
7. ```unassigned_locations.csv```<br>
This file contains the list of spots (locations) where no cells have been assigned by the algorithm. The columns include the spot IDs (rownames), the row (`row`) and column (`col`) indices of the spots, and the total number of cells in each spot (`Number of cells`), which is 0 for all the spots in this file.
8. ```log.txt```<br>
This file contains a log of CytoSPACE run parameters and running time.
</details>

## Frequently asked questions (FAQ)

<details><summary>Expand section</summary>

1. My ST dataset comes from a platform other than 10x Visium. Which additional parameters should I specify?<br>
Please refer to further instructions in the corresponding sections below:<br>
[__Running CytoSPACE on legacy ST data__](#running-cytospace-on-legacy-st-data) and [__Running CytoSPACE on single-cell ST data__](#running-cytospace-on-single-cell-st-data).

2. My scRNA-seq dataset comes in a format other than a UMI count matrix.<br>
You may want to make the following two changes to the default CytoSPACE workflow.<br>
(1) CytoSPACE's internal algorithm for estimating cell type fractions was generally written with a UMI count matrix in mind. As an alternative to using the internal algorithm, you may provide a `--cell-type-fraction-estimation-path` file instead. Please see [__Advanced Options__](#advanced-options) - __User-provided fractional composition of each cell type__ for more information.<br>
(2) CytoSPACE (v1.0.4+) downsamples scRNA-seq datasets to a certain number of transcripts per cell (1500 by default) prior to assignment so that the assignment is not dependent on the total transcript count of each cell. You can turn this feature off by appending a `--downsample-off` flag to the CytoSPACE call, and may instead choose to provide previously downsampled scRNA-seq data as input.

3. During runtime, I get a `Killed` / `Terminated` / `Segmentation fault` / `concurrent.futures.process.BrokenProcessPool` error.<br>
While these errors could arise from a variety of reasons, it is likely that the CytoSPACE run needs more memory than what is available.<br>
We provide a subsampling routine where the ST datasets are partitioned into smaller subsets and evaluated one subset at a time, which reduces memory requirements. Please see [__Advanced Options__](#advanced-options) - __Spot subsampling for parallelization__ for more information.<br>
If you are experiencing this error with a single-cell ST dataset, it will be helpful to reduce the `-noss` parameter instead.

4. My input data are very sparse. Is there a way to provide sparse matrices as input?<br>
As of CytoSPACE v1.1.0, CytoSPACE supports sparse matrices as input, for both single cell and spatial datasets. This applies both to the main cytospace function and to the get_cellfracs_seuratv3 R script for estimation of the fractional composition of each cell type. While the input still needs to follow a fixed logic of file naming, we broadened compatibility to include the most common naming conventions and extensions for different versions of Cell Ranger and Space Ranger sparse matrix outputs. Please refer to the [__Input Files__](#input-files) - __Using sparse matrices for gene expression__ section for details.
<!-- As of CytoSPACE v1.0.5, CytoSPACE supports sparse matrices as input. This requires a specific set of files to represent the sparse matrix, and we recommend that the R helper scripts under `Prepare_input_files` to be used for generating these input. Please refer to the [__Input Files__](#input-files) - __Using sparse matrices for gene expression__ section for details. -->

<!--
4. Is there a gene expression matrix for the results?<br>
We currently do not provide the gene expression matrix as an output file, as it is often very large in size and takes a long time to write to disk. However, the `OriginalCID` column of the `assigned_locations.csv` file will consist of the single-cell IDs from the input scRNA-seq expression matrix, which can be used to recover the gene expressions for each cell assigned to each spot.
-->

5. Providing the output from Space Ranger (v2.0.0+) results in an error.<br>
As of CytoSPACE v1.1.0, CytoSPACE supports the output from Space Ranger (v2.0.0+). Please refer to the [__Input Files__](#input-files) - __From Space Ranger outputs__ section for details.

<!-- We were notified that the instructions for providing the Space Ranger outputs directly as a tarball resulted in errors for the newer versions of Space Ranger. It seems that this is occuring due to a recent format change in Space Ranger outputs, and we are currently working to fix this issue. In the meantime, please use the standard four-file input format, with the ST gene expression and coordinates provided as two separate .txt files. We appreciate the users letting us know of the issue! -->

</details>

## Example dataset for running CytoSPACE

<details><summary>Expand section</summary>

For users to test CytoSPACE, we have included files for an example run:
- A HER2+ breast cancer scRNA-seq atlas by Wu et al. (<a href="https://www.nature.com/articles/s41588-021-00911-1" target="_blank">Nature Genetics, 2021</a>) and a HER2+ breast cancer FFPE specimen profiled by the Visium platform (<a href="https://www.10xgenomics.com/resources/datasets/human-breast-cancer-ductal-carcinoma-in-situ-invasive-carcinoma-ffpe-1-standard-1-3-0" target="_blank">10x Genomics</a>). Default parameters were selected with Visium samples in mind and are appropriate here.


### Download example datasets
A zip file containing the example dataset can be downloaded from the following link:
- <a href="https://drive.google.com/file/d/1G8gK4MxCmRG4JZi588wloMsP8iZlQf_z/view?usp=share_link" target="_blank">Breast cancer</a>
   
### Command for running example analysis:
Once the example files are downloaded and unzipped, the commands below can be run from inside the unzipped directory:
```bash
  cytospace -sp brca_scRNA_GEP.txt -ctp brca_scRNA_celllabels.txt -stp brca_STdata_GEP.txt -cp brca_STdata_coordinates.txt -o cytospace_results_brca -sm lap_CSPR
```
Please note that here we use the `lap_CSPR` solver for compatibility. If your system supports AVX2 intrinsics, you can run the same commands without the final argument to use the `lapjv` solver instead. __The CytoSPACE run should take around 5 minutes.__

### CytoSPACE output files for example breast cancer data
The main output from a CytoSPACE run is the file named `assigned_locations.csv`, which provides the ST spots to which the single cells have been assigned. 

<p align="center">
  <img width="600" src="images/assigned_locations.png">
</p>

The CytoSPACE results are visualized in heatmaps saved as `cell_type_assignments_by_spot.pdf` showing the distribution of single cells across ST spots for each cell type. Color bars indicate the number of cells of the respective cell type inferred per spot. Below are the heatmaps produced for the example BRCA data.

<p align="center">
  <img width="600" src="images/BRCA_cell_type_assignments_by_spot.png">
</p>

For comparison, consider the pathologist annotations of this ST sample as provided by 10x:

<p align="center">
  <img width="600" src="images/Visium_FFPE_Human_Breast_Cancer_Pathologist_Annotations.png">
</p>

CytoSPACE also provides a scatterplot showing cells of all types at once near their spot location, saved as `cell_type_assignments_by_spot_jitter.pdf`. Each cell is colored by its cell type. Below is the plot produced for the example BRCA data.

<p align="center">
  <img width="600" src="images/BRCA_cell_type_assignments_by_spot_jitter.png">
</p>

The number of cells per spot by cell type as well as in total are provided in the file `cell_type_assignments_by_spot.csv`. Fractional abundances of each cell type are returned in the file `fractional_abundances_by_spot.csv`. A log file recording CytoSPACE inputs and running times is output in the file `log.txt`.

A zip file of the expected CytoSPACE outputs (with CytoSPACE v1.0.0 using the `lap_CSPR` solver) are available to download at the following link:
- <a href="https://drive.google.com/file/d/1ZMA0XEl_pjC12mb8bZL8zI9yzYd2djbq/view?usp=share_link" target="_blank">Breast cancer results</a>

<details><summary><b>Simulated datasets</b></summary>

In addition to the example dataset above, the simulated datasets that we have generated for evaluating robustness of CytoSPACE across different conditions are available for download below.
They were generated using annotated Slide-seq datasets of mouse cerebellum and hippocampus sections from Rodriques et al. (<a href="https://www.science.org/doi/10.1126/science.aaw1219" target="_blank">Science, 2019</a>). Each simulated dataset contains subdirectories for data generated using different spot resolutions (5, 15, and 30 cells per spot), as well as an `scRNA` subdirectory containing reference single-cell datasets with perturbations in a defined percentage of genes. For more information, please see the Methods section of the [paper](https://www.nature.com/articles/s41587-023-01697-9).
1. <a href="https://drive.google.com/file/d/1qfz2T8u3HRG4qdZc9qafcO4aCvjA91Rb/view?usp=share_link" target="_blank">Cerebellum</a>
2. <a href="https://drive.google.com/file/d/1Jyd14n-ISc5lF65pnJWLhCCgSkpjtbsr/view?usp=share_link" target="_blank">Hippocampus</a>

</details>

</details>

## Running CytoSPACE on legacy ST data

<details><summary>Expand section</summary>

By default, the CytoSPACE parameters have been optimized for standard 10x Visium spatial slides. Datasets generated by the legacy ST platform can be run with similar commands, but we recommend that the following parameters be adjusted:
1. `--mean_cell_numbers`, or `-mcn`, should be set to `20`. The legacy ST platform has larger spot sizes, so we recommend mapping an average of 20 cells per spot.
2. `--geometry`, or `-g` should be set to `square`. This will allow the plot function to shape each spot as a square rather than a hexagon.

Similar to the example breast cancer dataset above, we provide an example dataset below:
- A melanoma scRNA-seq atlas by Tirosh et al. (<a href="https://www.science.org/doi/10.1126/science.aad0501?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed" target="_blank">Science, 2016</a>), and a melanoma specimen profiled by the legacy ST platform (Thrane et al., <a href="https://aacrjournals.org/cancerres/article/78/20/5970/631815/Spatially-Resolved-Transcriptomics-Enables" target="_blank">Cancer Research, 2018</a>).

The zip file containing the dataset can be downloaded <a href="https://drive.google.com/file/d/1hwK_sh355chdmW50yrPJq7_W8j6HuRHh/view?usp=share_link" target="_blank">here</a>.

Running CytoSPACE with the command below generates the results shown <a href="https://drive.google.com/file/d/1bX4SqrYzIXov_A5ivlJ8U0qD8_lXmmBf/view?usp=share_link" target="_blank">here</a> . The format of the output will be the same as the breast cancer dataset above. Please note that here we specify the `-ctfep` parameter instead of using CytoSPACE's internal algorithm for estimating cell fractions (see [__Advanced options__](#advanced-options) - __User-provided fractional composition of each cell type__) as the scRNA-seq atlas used as reference was generated using Smart-seq2. If using CytoSPACE v1.0.4+, the `--downsample-off` flag should additionally be specified.
```bash
  cytospace -sp melanoma_scRNA_GEP.txt -ctp melanoma_scRNA_celllabels.txt -stp melanoma_STdata_slide1_GEP.txt -cp melanoma_STdata_slide1_coordinates.txt -ctfep melanoma_cell_fraction_estimates.txt -o cytospace_results_melanoma -mcn 20 -g square -sm lap_CSPR
```
</details>

## Running CytoSPACE on single-cell ST data

<details><summary>Expand section</summary>

While designed for Visium-type data in which most spots contain RNA from multiple cells, CytoSPACE can also be used with single-cell resolution spatial data such as <a href="https://vizgen.com/resources/meet-the-merscope-platform/" target="_blank">Vizgen's MERSCOPE platform</a>. We expect this extension to be useful for reducing noise and expanding transcriptome coverage of each cell in the ST data, which in turn could allow for identifying spatially-dependent changes across genes more diverse than what a typical single-cell resolution ST platform alone can provide.

For the single-cell resolution mode, CytoSPACE partitions the ST data into smaller subsets and utilizes multiple CPU cores to assign downsampled versions of the reference scRNA-seq data to these regions.

To run CytoSPACE in single-cell mode, please add the following parameters to your command:
1. `--single-cell` (`-sc`), which tells CytoSPACE to run in single-cell mode.
2. Cell types for the ST dataset. Please note that for the single-cell mode, CytoSPACE does not support the internal estimation of cell type fraction, and the users are expected to specify one of the two options below.<br>
(1) `--st-cell-type-path` (`-stctp`)<br>
__If the cell types for individual spots are available, we recommend using this option.__ This file will list the cell type labels for each spot, in the same format as the scRNA-seq cell type labels specified under `--cell-type-path`. All of the cell types present in `--st-cell-type-path` must also be present in `--cell-type-path`.<br>
(2) `--cell-type-fraction-estimation-path` (`-ctfep`)<br>
However, if the user does not have access to the cell types for each individual spot, they can instead use this option. See the [__Advaced Options__](#advanced-options) - __User-provided fractional composition of each cell type__ section regarding how this file should be formatted. 
3. `--number-of-processors` (`-nop`), which denotes the number of cores to use.
4. `--number-of-selected-spots` (`-noss`), which denotes the number of ST spots in each subset. We generally recommend `-noss 10000`.

To run CytoSPACE with single-cell resolution spatial data:
 ```bash
 cytospace --single-cell \
    --scRNA-path /path/to/scRNA_geneexpression \
    --cell-type-path /path/to/scRNA_celllabels \
    --st-path /path/to/ST_geneexpression \
    --coordinates-path /path/to/ST_coordinates \
    --st-cell-type-path /path/to/ST_celllabels \
    --number-of-processors NUMBER_OF_PROCESSORS \
    --number-of-selected-spots NUMBER_OF_SELECTED_SPOTS
```
Or with more condensed parameter names: 
 ```bash
 cytospace -sc \
    -sp /path/to/scRNA_geneexpression \
    -ctp /path/to/scRNA_celllabels \
    -stp /path/to/ST_geneexpression \
    -cp /path/to/ST_coordinates \
    -stctp /path/to/ST_celllabels \
    -nop NUMBER_OF_PROCESSORS \
    -noss NUMBER_OF_SELECTED_SPOTS
```

A zip file of example single cell inputs is available to download from Google Drive <a href="https://drive.google.com/file/d/1odOcIfY3oqvLCNdXHLRaSmTraRxqnHLp/view?usp=share_link" target="_blank">here</a>.

To run CytoSPACE with this example dataset, run the following command from the location of the unzipped inputs and with your CytoSPACE conda environment active:
 ```bash
 cytospace \
    -sp HumanColonCancerPatient2_scRNA_expressions_cytospace.tsv \
    -ctp HumanColonCancerPatient2_scRNA_annotations_cytospace.tsv \
    -stp HumanColonCancerPatient2_ST_expressions_cytospace.tsv \
    -cp HumanColonCancerPatient2_ST_coordinates_cytospace.tsv \
    -stctp HumanColonCancerPatient2_ST_celltypes_cytospace.tsv \
    -o cytospace_results_crc \
    -sm lap_CSPR \
    -sc -noss 10000 -nop 2
```

Running CytoSPACE in the `--single-cell` mode will output the assignments `assigned_locations.csv`, the plot `cell_type_assignments_by_spot_single_cell.pdf`, and the log file `log.txt`. The plot generated will be a scatterplot of the cells colored by cell type, as shown below for the example dataset. The full results for the example dataset using the above command is available for download <a href="https://drive.google.com/file/d/1LTTDVGAuQ4QYkyCX6WtyBNXcnZe9fxKG/view?usp=share_link" target="_blank">here</a>.

<p align="center">
  <img width="600" src="images/CRC_cell_type_assignments_by_spot_single_cell.png">
</p>

</details>

## Advanced options
While the default options are recommended for most use cases, we do provide additional advanced options.

<details><summary><b>User-provided fractional composition of each cell type (-ctfep)</b></summary>

To account for the disparity between scRNA-seq and ST data in the number of cells per cell type, CytoSPACE requires the fractional composition of each cell type in the ST tissue. By default, CytoSPACE will generate this information by internally calling the `get_cellfracs_seuratv3.R` script using the input files. This script uses `Seurat v3`, which is installed as part of the CytoSPACE environment. We highly recommend using `Seurat v3` over `Seurat v4` for the purposes of cell type fraction estimation.

While our provided script uses <a href="https://satijalab.org/seurat/articles/spatial_vignette.html" target="_blank">Spatial Seurat</a>, there is a diverse set of other approaches available such as <a href="https://www.sanger.ac.uk/tool/cell2location/" target="_blank">cell2location</a>, <a href="https://github.com/MarcElosua/SPOTlight" target="_blank">SPOTlight</a>, or <a href="https://cibersortx.stanford.edu/" target="_blank">CIBERSORTx</a>.

Users may choose to skip CytoSPACE's internal algorithm and instead provide their own file for the estimated cell type composition of the ST dataset, specified with the `--cell-type-fraction-estimation-path` (`-ctfep`) flag. In particular, we recommend that a separate `-ctfep` file be provided if the reference scRNA-seq dataset comes from technologies that are not based on UMI counts, such as Smart-seq.

The provided file must be a table consisting of 2 rows with row names, where the first row contains the cell types and the second row contains the cell fractions of each cell type represented as proportions between 0 and 1. __Please make sure that the cell type labels in the first row match the labels present in the cell type label file, and that the cell type fractions sum to one. Row names must be present for both rows.__
<p align="center">
  <img src="images/cell_type_fractions_file.png">
</p>
</details>

<details><summary><b>User-provided estimates of number of cells per spot (-ncpsp)</b></summary>

Rather than using the internal mechanism of CytoSPACE for estimating the number of cells per spot, users can provide their own estimates (from image segmentation, for example) in a two-column file with header, in which the first column contains spot IDs and the second contains the number of cells predicted per spot:

<p align="center">
  <img width="300" src="images/n_cells_per_spot.PNG">
</p>

To run CytoSPACE with this option, pass the flag `-ncpsp` or `--n-cells-per-spot-path` followed by the file location.
</details>

<details><summary><b>Spot subsampling for parallelization (-sss)</b></summary>

The memory and runtime required for running CytoSPACE may vary based on the number of spots. To allow for CytoSPACE to run under different conditions, we provide an option to partition the estimated number of cells in the ST sample into smaller chunks, where similarly downsampled reference scRNA-seq data are then assigned using multiple CPU cores. Please note that this option is only for non-single-cell ST datasets (10x Visium and legacy ST); users running CytoSPACE on single-cell ST data should instead modify the `-noss` parameter to achieve the same effect.

The users can use this option by specifying the `--sampling-sub-spots` (`-sss`) flag, along with the desired number of subsampled cells per partition (`--number-of-selected-sub-spots`, or `-nosss`) and the number of cores to be used (`--number-of-processors`, or `-nop`). Please note that the right `-nosss`/`-nop` combination will be highly dependent on the user's system configurations, and one may need to try different combinations to see which one allows for a successful run. Ideally, using the highest possible `-nosss` without going out of memory would be the most effective. We noticed that `-sss -nosss 3000 -nop 2` works well for the example BRCA dataset on an environment with 8GB RAM.

The following example command will run CytoSPACE on the example breast cancer dataset, assigning scRNA-seq data to 5000 cells at a time using 2 cores:
```bash
  cytospace \
    -sp brca_scRNA_GEP.txt \
    -ctp brca_scRNA_celllabels.txt \
    -stp brca_STdata_GEP.txt \
    -cp brca_STdata_coordinates.txt \
    -o cytospace_results_brca \
    -sm lap_CSPR \
    -sss -nosss 5000 -nop 2
```
</details>

<details><summary><b>Alternative distance metric (-dm)</b></summary>

By default, CytoSPACE uses Pearson correlation to compare cell and spot transcriptomes. Users can choose to use Spearman correlation or Euclidean distance instead by passing `-dm Spearman_correlation` or `-dm Euclidean` respectively with the function call. 
</details>

<details><summary><b>Setting a new random seed (-se)</b></summary>

While the CytoSPACE algorithm is mostly deterministic, the initial step of sampling cells to be mapped is done at random. To provide an alternative random seed resulting in a different random sampling of cells, users can pass `-se` followed by the desired (integer) seed with the function call. The default random seed for CytoSPACE is 1.
</details>

<details><summary><b>Alternative handling of sampling (-sam)</b></summary>

CytoSPACE starts by creating a pool of cells that matches what is expected within the ST data. By default, this is done by resampling single cells to achieve the overall cell type fractions and total cell numbers estimated in the tissue. We recommend that CytoSPACE be run with this default setting for all real data analyses. However, we provide an additional option to generate new "place-holder" cells by sampling from the distribution of gene counts within each cell type instead, and used this option for ensuring uniqueness of mapped cells for benchmarking on simulated data. To run CytoSPACE with this alternative mode, users can pass `-sam place_holders` with the function call. When running in place-holder mode, the `assigned_expression` directory will not be generated; instead, the gene expression of the newly generated "place-holder" cells will be saved as part of the output under `new_scRNA.csv`.
</details>

<details><summary><b>Method extension: mapping quality</b></summary>

While CytoSPACE's formulation as a linear assignment problem guarantees an optimal solution given its cost function, there is no underlying probabilistic framework for estimating mapping uncertainty. One possibility is to determine whether a given cell type belongs to a given spot after mapping - that is, whether a spot contains at least one cell of the same cell type. Notably, this does not distinguish between cells of the same cell type for quality of fit. As such a protocol provides some measure of mapping quality, albeit incomplete, we provide a helper script that implements this via a support vector machine that produces and trains on pseudo-bulks generated from the input scRNA-seq data. This script, `uncertainty_quantification.R`, takes as input the path to the ST dataset count matrix file, the scRNA-seq count matrix file, and the CytoSPACE output file `assigned_locations.csv`, and returns an appended output file with confidence scores in `assigned_locationswConfidenceScores.csv`. The command to run this script following a completed CytoSPACE run is as follows: 
 ```bash
 Rscript uncertainty_quantification.R /path/to/ST_geneexpression /path/to/scRNA_geneexpression /path/to/assigned_locations.csv
```
For interpreting confidence scores, we recommend a cutoff of 0.1, with higher scores indicating increased confidence that a spot contains at least one cell of the same cell type.

Please note that `uncertainty_quantification.R` requires separate dependencies from those included in the provided `environment.yml` file for the `cytospace` conda environment. This script should be run in a separate environment with the following R packages installed: `Seurat` (must be v4; tested with v4.0.1), `data.table` (tested with v1.14.0), and `e1071` (tested with v1.7.8).
</details>

## Extended usage details

<details><summary>Expand section</summary>

```
usage: cytospace [-h] -sp SCRNA_PATH -ctp CELL_TYPE_PATH [-stp ST_PATH] [-cp COORDINATES_PATH] [-srp SPACERANGER_PATH]
                 [-stctp ST_CELL_TYPE_PATH] [-ctfep CELL_TYPE_FRACTION_ESTIMATION_PATH] [-ncpsp N_CELLS_PER_SPOT_PATH]
                 [-o OUTPUT_FOLDER] [-op OUTPUT_PREFIX] [-mcn MEAN_CELL_NUMBERS] [--downsample-off]
                 [-smtpc SCRNA_MAX_TRANSCRIPTS_PER_CELL] [-sc] [-noss NUMBER_OF_SELECTED_SPOTS] [-sss]
                 [-nosss NUMBER_OF_SELECTED_SUB_SPOTS] [-nop NUMBER_OF_PROCESSORS] [-sm {lapjv,lapjv_compat,lap_CSPR}]
                 [-dm {Pearson_correlation,Spearman_correlation,Euclidean}] [-sam {duplicates,place_holders}]
                 [-se SEED] [-p] [-g GEOMETRY] [-nc NUM_COLUMN] [-mp MAX_NUM_CELLS_PLOT]

CytoSPACE is a computational strategy for assigning single-cell transcriptomes to in situ spatial transcriptomics (ST)
data. Our method solves single cell/spot assignment by minimizing a correlation-based cost function through a linear
programming-based optimization routine.

optional arguments:
  -h, --help            show this help message and exit
  -stp ST_PATH, --st-path ST_PATH
                        Path to spatial transcriptomics data (expressions)
  -cp COORDINATES_PATH, --coordinates-path COORDINATES_PATH
                        Path to transcriptomics data (coordinates)
  -srp SPACERANGER_PATH, --spaceranger-path SPACERANGER_PATH
                        Path to SpaceRanger tar.gz data file
  -stctp ST_CELL_TYPE_PATH, --st-cell-type-path ST_CELL_TYPE_PATH
                        Path to ST cell type file (recommended for single-cell ST)
  -ctfep CELL_TYPE_FRACTION_ESTIMATION_PATH, --cell-type-fraction-estimation-path CELL_TYPE_FRACTION_ESTIMATION_PATH
                        Path to ST cell type fraction file (recommended for bulk ST)
  -ncpsp N_CELLS_PER_SPOT_PATH, --n-cells-per-spot-path N_CELLS_PER_SPOT_PATH
                        Path to number of cells per ST spot file
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        Relative path to the output folder
  -op OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Prefix of results stored in the 'output_folder'
  -mcn MEAN_CELL_NUMBERS, --mean-cell-numbers MEAN_CELL_NUMBERS
                        Mean number of cells per spot, default 5 (appropriate for Visium). If analyzing legacy spatial
                        transcriptomics data, set to 20
  --downsample-off      Turn off downsampling for scRNA-seq data
  -smtpc SCRNA_MAX_TRANSCRIPTS_PER_CELL, --scRNA_max_transcripts_per_cell SCRNA_MAX_TRANSCRIPTS_PER_CELL
                        Number of transcripts per cell to downsample scRNA-seq dataset to. This allows for assignments
                        that are not dependent on the overall expression level
  -sc, --single-cell    Use single-cell spatial approach if specified
  -noss NUMBER_OF_SELECTED_SPOTS, --number-of-selected-spots NUMBER_OF_SELECTED_SPOTS
                        Number of selected spots from ST data used in each iteration
  -sss, --sampling-sub-spots
                        Sample subspots to limit the number of mapped cells if specified
  -nosss NUMBER_OF_SELECTED_SUB_SPOTS, --number-of-selected-sub-spots NUMBER_OF_SELECTED_SUB_SPOTS
                        Number of selected subspots from ST data to limit the number of mapped cells
  -nop NUMBER_OF_PROCESSORS, --number-of-processors NUMBER_OF_PROCESSORS
                        Number of processors used for the analysis
  -sm {lapjv,lapjv_compat,lap_CSPR}, --solver-method {lapjv,lapjv_compat,lap_CSPR}
                        Which solver to use for the linear assignment problem, default 'lapjv'
  -dm {Pearson_correlation,Spearman_correlation,Euclidean}, --distance-metric {Pearson_correlation,Spearman_correlation,Euclidean}
                        Which distance metric to use for the cost matrix, default 'Pearson_correlation'
  -sam {duplicates,place_holders}, --sampling-method {duplicates,place_holders}
                        Which underlying method to use for dealing with duplicated cells, default 'duplicates'
  -se SEED, --seed SEED
                        Set seed for random generators, default 1
  -p, --plot-off        Turn create plots on/off
  -g GEOMETRY, --geometry GEOMETRY
                        ST geometry, either 'honeycomb' or 'square' accepted
  -nc NUM_COLUMN, --num-column NUM_COLUMN
                        Number of columns in figure
  -mp MAX_NUM_CELLS_PLOT, --max-num-cells-plot MAX_NUM_CELLS_PLOT
                        Maximum number of cells to plot in single-cell visualization

Required arguments:
  -sp SCRNA_PATH, --scRNA-path SCRNA_PATH
                        Path to scRNA-Seq data
  -ctp CELL_TYPE_PATH, --cell-type-path CELL_TYPE_PATH
                        Path to cell type labels
```

You can see this list of variables and default values for running CytoSPACE from the commmand line as well at any time by calling `cytospace` along with the `-h` or 
`--help` flag, i.e., `cytospace -h`.
</details>

## CytoSPACE solver options

<details><summary>Expand section</summary>

1. `lapjv` __(Recommended for most systems)__    By default, CytoSPACE calls the `lapjv` solver from package `lapjv`. This solver is a fast implementation of the Jonker-Volgenant shortest augmenting path assignment algorithm and returns a globally optimal solution given the objective function as defined in our [paper](https://www.nature.com/articles/s41587-023-01697-9). As noted above, however, this package is not supported on all systems as it achieves its speedup through use of AVX2 instructions. This solver will be selected by default and can be specified explicitly by passing arguments `--solver-method lapjv` or `-sm lapjv` to `cytospace`.
2. `lap_CSPR` __(Recommended for systems not supporting `lapjv`)__    A second solver option is the `linear_assignment` method from the `ortools` package. This solver uses a different method than the first and third options, an assignment algorithm called the cost scaling push relabel method. This algorithm approximates assignment costs to integer values and loses some numerical precision in doing so. Therefore, while it returns a globally optimal solution __after approximation__ given the objective function defined in the paper, it will return similar but generally not identical results to the first two methods. This solver has a similar running time to the first option and is a good option for systems not supporting the `lapjv` package. This solver can be selected by passing arguments `--solver-method lap_CSPR` or `-sm lap_CSPR` to `cytospace`.
3. `lapjv_compat`   A third solver option implements the `lapjv` solver from package `lap`. Like the first option `lapjv`, this solver also implements the Jonker-Volgenant shortest augmenting path assignment algorithm to return the same globally optimal solution given the objective function defined in the paper. Furthermore, it is broadly supported and should work on all standard operating systems. However, it takes 3-4 times as long to run as the first solver option (the `lapjv` solver from the `lapjv` package), so we only recommend it for systems that do not support the first option. This solver can be selected by passing arguments `--solver-method lapjv_compat` or `-sm lapjv_compat` to `cytospace`.
</details>

## Updating local installations

<details><summary>Expand section</summary>

To update your local installation of CytoSPACE following updates of this GitHub repository, navigate to your `cytospace` directory and execute the following commands:
```bash
  git pull
  conda env update --name cytospace --file environment.yml
  conda activate cytospace
  pip install .
```
If you have made local updates to your version of the CytoSPACE source code, you should execute 
```bash
  pip install .
``` 
before running. 
</details>

## Authors
CytoSPACE was developed in the <a href="https://anlab.stanford.edu/" target="_blank">Newman Lab</a> by

* Milad R. Vahid (miladrv)
* Erin L. Brown (erinlbrown)
* Chloé B. Steen (cbsteen)
* Wubing Zhang (WubingZhang)
* Hyun Soo Jeon (hsjeon-k)
* Aaron M. Newman (aaronmnewman)

## Contact
If you have any questions, please contact the CytoSPACE team at cytospaceteam@gmail.com.

## License
Please see the <a href="LICENSE" target="_blank">LICENSE</a> file.

## Citation
If you use CytoSPACE, please cite:  

*High-resolution alignment of single-cell and spatial transcriptomes with CytoSPACE* (Nature Biotechnology 2023)  
Milad R. Vahid*, Erin L. Brown*, Chloé B. Steen*, Wubing Zhang, Hyun Soo Jeon, Minji Kang, Andrew J. Gentles, Aaron M. Newman.  
https://www.nature.com/articles/s41587-023-01697-9
