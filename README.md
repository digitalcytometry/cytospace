<p align="center">
  <img width="300" src="https://github.com/digitalcytometry/cytospace/blob/main/images/CytoSPACE_logo.jpeg">
</p>

<h1> <p align="center">
    Robust and rapid alignment of single-cell and spatial transcriptomes
</p> </h1>

**CytoSPACE** is a novel computational tool for assigning single-cell transcriptomes to in situ spatial transcriptomics (ST) data. Our method solves single cell/spot assignment by minimizing a correlation-based cost function through a shortest augmenting path optimization routine. 

<p align="center">
  <img src="https://github.com/digitalcytometry/cytospace/blob/main/images/CytoSPACE_overview.png" width="900"> 
</p>

The key innovations of our method are:

- Unlike conventional methods which calculate cell type decompositions by spot, CytoSPACE yields a reconstructed tissue specimen with both high gene coverage and spatially-resolved scRNA-seq data suitable for downstream analysis.
- CytoSPACE is highly robust to noise and returns globally optimal cell-to-spot assignments. (See the paper for full details.)
- Unlike other methods which generally operate on pre-selected marker genes or on a shared embedding space (the latter of which can erase true biological variation), CytoSPACE uses the full transcriptome without the need for batch correction, helping it retain sensitivity to subtle cell states.
- CytoSPACE is quick and simple to execute. It runs in minutes even with a single CPU on a personal laptop and requires no hyperparameter tuning or gene/feature selection.

Overview of README:
-   [**Installation instructions**](#installation-instructions)
-   [**File formatting**](#file-formatting)
-   [**Preprocessing: Input file preparation**](#preprocessing-input-file-preparation)
-   [**Preprocessing: Cell fraction estimation with Seurat v3**](#preprocessing-cell-fraction-estimation-with-seurat-v3)
-   [**Running CytoSPACE**](#running-cytospace)
-   [**CytoSPACE outputs**](#cytospace-outputs)
-   [**Example datasets for running CytoSPACE**](#example-datasets-for-running-cytospace)
-   [**Customizing plotting outputs**](#customizing-plotting-outputs)
-   [**Advanced options**](#advanced-options)
-   [**Extended usage details**](#extended-usage-details)
-   [**CytoSPACE solver options**](#cytospace-solver-options)
-   [**Updating local installations**](#updating-local-installations)
-   [**Authors**](#authors)
-   [**Contact**](#contact)
-   [**License**](#license)
-   [**Citation**](#citation)


## Installation instructions (5-10 minutes)
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
We highly recommend you install this package, which provides a fast implementation of the default core optimization algorithm within CytoSPACE. However, some systems may not accommodate it as it requires CPU support for AVX2 instructions. To determine if your system supports this package, it is generally easiest to simply attempt to install it as above. If it installs without problems, your system will support it! If you run into an error, it is likely your system does not support it, and you can simply use one of the other options we have provided. See __Solver options__ below for details. Please note that if the package installs but you receive an "illegal instruction" error while running CytoSPACE, you may be able to build the package instead with the following command:
```bash
   pip3 install git+https://github.com/src-d/lapjv
```
For more information, see the <a href="https://pypi.org/project/lapjv/" target="_blank">lapjv documentation page</a>. 

## File formatting
CytoSPACE requires 5 files as input. All files should be provided in tab-delimited tabular input format (saved as .txt) with no double quotations. Further formatting details for each input file are specified below. For helper scripts to prepare input files directly from Seurat objects, refer to the section [Preprocessing: Input file preparation](#preprocessing-input-file-preparation).

1. __A scRNA-seq gene expression file:__
- The matrix must be genes (rows) by cells (columns).
- The first row must contain the single cell IDs and the first column must contain the gene names.
- The first column (gene names) must have a header.
- The gene expression data should be represented as non-normalized counts. 
<img src="https://github.com/digitalcytometry/cytospace/blob/main/images/scRNAfile.png" width="800"> 

2. __A cell type label file:__
- Cell type labels corresponding to the single cell IDs in the scRNA-seq gene expression matrix. 
- Cell type label strings should not have special characters. 
- The table should contain two columns, where column 1 contains the single cell IDs corresponding to the columns of the scRNA-seq matrix and column 2 contains the corresponding cell type labels.
- The columns must have a header. 
<img src="https://github.com/digitalcytometry/cytospace/blob/main/images/celllabelfile.png" width="250"> 

3. __A spatial transcriptomics (ST) gene expression file:__
- The matrix must be genes (rows) by ST spots (columns).
- The first row must contain the ST spot IDs and the first column must contain the gene names.
- The first column (gene names) must have a header.
- The gene expression data should be represented as non-normalized counts. 
<img src="https://github.com/digitalcytometry/cytospace/blob/main/images/STdatafile.png" width="800"> 

4. __A spatial transcriptomics coordinates file:__
- A table consisting of 3 columns, where the first column contains the ST spot IDs corresponding to and in the same order as the columns of the ST gene expression matrix, and column 2 and 3 contain the row and column indices of the spatial transcriptomics data, respectively. 
- The columns must have a header. 
<img src="https://github.com/digitalcytometry/cytospace/blob/main/images/STcoordfile.png" width="250"> 

5. __A file with cell type fraction estimates:__
- A table consisting of 2 rows, where the first row is the cell type labels, and the second row is the cell fractions of each cell type represented as proportions between 0 and 1. The first column is the row names.
- __The first row of cell type labels must match the labels present in the cell type label file.__
- The cell type fractions should sum to one.  
- Cell type fractions can be generated via any method for fractional abundance estimation. We provide one such implementation using Spatial Seurat in `get_cellfracs_seuratv3.R`. For further details, see section "__Preprocessing__" below.
<img src="https://github.com/digitalcytometry/cytospace/blob/main/images/cell_type_fractions_file.png" width="800"> 

                                                                                                                 
## Preprocessing: Input file preparation
If you have data in the form of Seurat objects, you can generate files formatted for CytoSPACE input via helper functions we have provided in the `R` script `generate_cytospace_from_seurat_object.R` in `cytospace/Prepare_input_files`. To use these helper functions, first import them from `generate_cytospace_from_seurat_object.R` by including 
```bash
  source('/path/to/generate_cytospace_from_seurat_object.R')
```
in your R script. 

### From scRNA-seq Seurat object
For producing CytoSPACE inputs from scRNA Seurat objects, we provide the function `generate_cytospace_from_scRNA_seurat_object` which may be called as
```bash
  generate_cytospace_from_scRNA_seurat_object(scRNA_Seurat_Object,dir_out='',fout_prefix='')
```
within your R script. The first argument (required) designates your input Seurat object, `dir_out` (optional, default is working directory) specifies the path to the output directory to store the results, and `fout_prefix` (optional, default is none) specifies a prefix to add to output file names, which otherwise are generated as `scRNA_data.txt` and `cell_type_labels.txt`. Please note that `Idents(scRNA_Seurat_Object)` must be set to include cell types.


### From Spatial Seurat object
For producing CytoSPACE inputs from ST Seurat objects, we provide the function `generate_cytospace_from_ST_seurat_object` which may be called as
```bash
  generate_cytospace_from_ST_seurat_object(ST_Seurat_Object,dir_out='',fout_prefix='',slice='slice1')
```
within your R script. The first argument (required) designates your input Seurat object, `dir_out` (optional, default is working directory) specifies the path to the output directory to store the results, `fout_prefix` (optional, default is none) specifies a prefix to add to output file names, which otherwise are generated as `ST_data.txt` and `Coordinates.txt`, and `slice` (optional, default is `slice1`) provides the name of your slice as stored in your Seurat object.

### From raw Space Ranger outputs
If you are starting from Space Ranger outputs, an easy way to prepare inputs formatted for CytoSPACE is by loading the Space Ranger outputs into a Seurat object with
<a href="https://satijalab.org/seurat/reference/load10x_spatial" target="_blank">Load10X_Spatial</a>. Once you have a Seurat object, you can use the helper script provided above.

## Preprocessing: Cell fraction estimation with Seurat v3
To account for the disparity between scRNA-seq and ST data in the number of cells per cell type, the fractional composition of each cell type in the ST tissue needs to be provided as input to CytoSPACE. This is determined using an external deconvolution tool, such as <a href="https://satijalab.org/seurat/articles/spatial_vignette.html" target="_blank">Spatial Seurat</a>, <a href="https://www.sanger.ac.uk/tool/cell2location/" target="_blank">cell2location</a>, <a href="https://github.com/MarcElosua/SPOTlight" target="_blank">SPOTlight</a>, or <a href="https://cibersortx.stanford.edu/" target="_blank">CIBERSORTx</a>. We have included Spatial Seurat in our benchmarking, and provide here a script to obtain the cell type fractions using this approach.

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
2. __You must run this script within the `cytospace` conda environment.__ This is done using the command `conda activate cytospace`.
3. We use `Seurat v3` for estimating cell fractions, and this is installed as part of the CytoSPACE environment. If you want to run other analyses using more recent versions of Seurat after running CytoSPACE, for example Seurvat v4, make sure to first __deactivate the CytoSPACE environment__ once you are done running CytoSPACE. This is done using the command `conda deactivate cytospace`.
4. __We highly recommend using Seurat v3 rather than v4 for cell type fraction estimation.__ Running cell type fraction estimation from within the `cytospace` environment will ensure this. If you would like to check your Seurat version within an R session, you can use the following command:
```
library(Seurat)
packageVersion('Seurat')
```

## Running CytoSPACE
After activating the `cytospace` conda environment via `conda activate cytospace`, CytoSPACE can be called from the command line from any folder using `cytospace`. Examples on how to run CytoSPACE are provided in the section "Example datasets for running CytoSPACE" below.

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

For full usage details with additional options, see the __Extended usage details__ section below. 

### Choosing a solver
CytoSPACE provides three solver options. In short, we recommend using the default option `lapjv` if your system supports AVX2 (i.e., if you were able to successfully install it with `pip install lapjv==1.3.14`) and `lap_CSPR` otherwise. No options are required to use the default solver `lapjv`. To use `lap_CSPR` instead, pass the argument `-sm lap_CSPR` to your `cytospace` call. For full solver details, see the __Solver options__ section below.


### Other ways CytoSPACE can be run:
1. You can call the `cytospace.py` script directly with python:
 `  python cytospace/cytospace.py`
 
2. You can import methods or functions from `CytoSPACE` in python and modify/create your own 
    pipeline. For example:
```python
  from cytospace import cytospace

  for mean_cell_numbers in [5, 10, 20]:
      cytospace(..., mean_cell_numbers=mean_cell_numbers)
```

## CytoSPACE outputs
CytoSPACE will produce five output files by default.
1. ```plot_cell_type_locations.pdf``` Heatmaps of cell type assignments within the ST sample. Along with a plot showing the total number of cells mapped to each spot, these show the spatial distribution of cell type assignments. Color bars indicate the number of cells of the respective cell type inferred per spot.
2. ```assigned_locations.csv``` This file will provide the assigned locations of each single cell mapped to ST spots. As some cells may be mapped to multiple locations depending on the size of the input scRNA-seq set, new cell IDs (`UniqueCID`) are assigned to each cell and given in the first column. The second column includes original cell IDs (`OriginalCID`); the third includes assigned spot IDs (`SpotID`); the fourth and fifth columns respectively include  `row` and `column` indices of the corresponding spots; and then optionally, the sixth and seventh columns include `coord_x` and `coord_y` of the corresponding spots if these details were provided in the initial Coordinates file.
3. ```cell_type_assignments_by_spot.csv``` This file gives the raw number of cells of each cell type per spot by `SpotID` as well as the total number of cells assigned to that spot.
4. ```fractional_abundances_by_spot.csv``` This file gives the fractional abundance of cell types assigned to each spot by `SpotID`.
5. ```log.txt``` This file contains a log of CytoSPACE run parameters and running time.

## Example datasets for running CytoSPACE
For users to test CytoSPACE, we have included files for two example runs:
1. A HER2+ breast cancer scRNA-seq atlas by Wu et al. (<a href="https://www.nature.com/articles/s41588-021-00911-1" target="_blank">Nature Genetics, 2021</a>) and a HER2+ breast cancer FFPE specimen profiled by the Visium platform (<a href="https://www.10xgenomics.com/resources/datasets/human-breast-cancer-ductal-carcinoma-in-situ-invasive-carcinoma-ffpe-1-standard-1-3-0" target="_blank">10x Genomics</a>). Default parameters were selected with Visium samples in mind and are appropriate here.
2. A melanoma scRNA-seq atlas by Tirosh et al. (<a href="https://www.science.org/doi/10.1126/science.aad0501?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed" target="_blank">Science, 2016</a>), and a melanoma specimen profiled by the legacy ST platform (Thrane et al, <a href="https://aacrjournals.org/cancerres/article/78/20/5970/631815/Spatially-Resolved-Transcriptomics-Enables" target="_blank">Cancer Research, 2018</a>). As the legacy ST platform has larger spot sizes, we recommend mapping an average of 20 cells per spot, i.e. passing the argument `-mcn 20` to your `cytospace` call.

### Download example datasets
Zip files containing the example datasets can be downloaded from the following links:
1. <a href="https://drive.google.com/file/d/1ODN7Ux2x9XNl1K8cOLl4JxIZpjJL2GwG/view?usp=sharing" target="_blank">Breast cancer</a>
2. <a href="https://drive.google.com/file/d/1oHe4UP2K0kQS9gNFVvZtjqtJyJK_VeBp/view?usp=sharing" target="_blank">Melanoma</a> 

To download from the command line using `gdown`:
1. Breast cancer
   ```bash
   gdown --fuzzy https://drive.google.com/file/d/1ODN7Ux2x9XNl1K8cOLl4JxIZpjJL2GwG/view?usp=sharing
   unzip CytoSPACE_example_breast_cancer.zip
   ```
2. Melanoma
   ```bash
   gdown --fuzzy https://drive.google.com/file/d/1oHe4UP2K0kQS9gNFVvZtjqtJyJK_VeBp/view?usp=sharing
   unzip CytoSPACE_example_melanoma.zip
   ```
   
### Commands for running example analyses:
Once the example files are downloaded, the commands below can be run from the folders where the example datasets are located:
```bash
  cytospace -sp brca_scRNA_GEP.txt -ctp brca_scRNA_celllabels.txt -stp brca_STdata_GEP.txt -cp brca_STdata_coordinates.txt -ctfep brca_cell_fraction_estimates.txt -o cytospace_results_brca -sm lap_CSPR
```

```bash
  cytospace -sp melanoma_scRNA_GEP.txt -ctp melanoma_scRNA_celllabels.txt -stp melanoma_STdata_slide1_GEP.txt -cp melanoma_STdata_slide1_coordinates.txt -ctfep melanoma_cell_fraction_estimates.txt -mcn 20 -o cytospace_results_melanoma -nr 5 -nc 3 -ss 1100 -pm s -nv -sm lap_CSPR 
```
Please note that here we use the `lap_CSPR` solver for compatibility. If your system supports AVX2 intrinsics, you can run the same commands without the final argument to use the `lapjv` solver instead. __These CytoSPACE runs should take around 5 minutes each.__

### CytoSPACE output files for example breast cancer data
The main output from a CytoSPACE run is the file named `assigned_locations.csv`, which provides the ST spots to which the single cells have been assigned. 

<img width="600" src="https://github.com/digitalcytometry/cytospace/blob/main/images/assigned_locations.png">


The CytoSPACE results are visualized in heatmaps saved as `plot_cell_type_locations.pdf` showing the distribution of single cells across ST spots for each cell type. Color bars indicate the number of cells of the respective cell type inferred per spot. Below are the heatmaps produced for the example BRCA data.

<p align="center">
  <img width="800" src="https://github.com/digitalcytometry/cytospace/blob/main/images/BRCA_plot_cell_type_locations.png">
</p>

For comparison, consider the pathologist annotations of this ST sample as provided by 10x:

<p align="center">
  <img width="800" src="https://github.com/digitalcytometry/cytospace/blob/main/images/Visium_FFPE_Human_Breast_Cancer_Pathologist_Annotations.png">
</p>

The number of cells per spot by cell type as well as in total are provided in the file `cell_type_assignments_by_spot.csv`. Fractional abundances of each cell type are returned in the file `fractional_abundances_by_spot.csv`. A log file recording CytoSPACE inputs and running times is output in the file `log.txt`.

Zip files of expected CytoSPACE outputs (with `lap_CSPR` solver) are available to download at the following links:
1. <a href="https://drive.google.com/file/d/1CLfy4Txez8ThID8YzIH04hlvrBRCQ4Rh/view?usp=sharing" target="_blank">Breast cancer results</a>
2. <a href="https://drive.google.com/file/d/1X4jMwctRNmqCRIJcop2hL3jRdlhxnnlc/view?usp=sharing" target="_blank">Melanoma results</a> 

To download from the command line using `gdown`:
1. Breast cancer
   ```bash
   gdown --fuzzy https://drive.google.com/file/d/1CLfy4Txez8ThID8YzIH04hlvrBRCQ4Rh/view?usp=sharing
   unzip CytoSPACE_example_breast_cancer_results.zip
   ```
2. Melanoma
   ```bash
   gdown --fuzzy https://drive.google.com/file/d/1X4jMwctRNmqCRIJcop2hL3jRdlhxnnlc/view?usp=sharing
   unzip CytoSPACE_example_melanoma_results.zip
   ```

## Customizing plotting outputs
CytoSPACE provides two mechanisms for plotting output heatmaps. First, CytoSPACE generates heatmaps by default within the main function call.  By default, the geometry and plotting parameters have been optimized for standard 10x Visium spatial slides. To plot non-Visium spatial data, the `plot-nonvisium` flag should be set to `True`. As an example, to change the plotting parameters, e.g., `-nr` (number of rows of heatmaps per page), `-nc` (number of columns of heatmaps per page), `-ss` (spot size), -pm (plot marker which indicates the shape of the spots, e.g., `-pm s` plots squares and `-pm h` plots hexagons) and `-nv` (plot non-Visium), from command line, for the melanoma dataset:
```bash
  cytospace -sp melanoma_scRNA_GEP.txt -ctp melanoma_scRNA_celllabels.txt -stp melanoma_STdata_slide1_GEP.txt -cp melanoma_STdata_slide1_coordinates.txt -ctfep melanoma_cell_fraction_estimates.txt -nr 5 -nc 3 -ss 1100 -pm s -nv
```
`cytospace-plot` can be used to generate heatmaps from CytoSPACE outputs. There are three paths that need to be set by users: `-alp` (assigned_locations_path), `-cp` (coordinates_path) and `-o` (output_filename). After CytoSPACE outputs have been generated, `cytospace-plot` can be run from the command line for the sample breast cancer data as:
```bash
  cytospace-plot -alp cytospace_results_brca/assigned_locations.csv -cp  brca_STdata_coordinates.txt -o brca_results.pdf
```
and for the sample melanoma data as: 
```bash
  cytospace-plot -alp cytospace_results_melanoma/assigned_locations.csv -cp  melanoma_STdata_slide1_coordinates.txt -o melanoma_results.pdf -nr 5 -nc 3 -ss 1100 -pm s -nv
```

## Advanced options
While default options are recommended for most use cases, we do provide additional advanced options.

### User-provided estimates of number of cells per spot
Rather than using the internal mechanism of CytoSPACE for estimating the number of cells per spot, users can provide their own estimates (from image segmentation, for example) in a two-column file with header, in which the first column contains spot IDs and the second contains the number of cells predicted per spot:

<p align="center">
  <img width="300" src="https://github.com/digitalcytometry/cytospace/blob/main/images/n_cells_per_spot.PNG">
</p>

To run CytoSPACE with this option, pass the flag `-ncpsp` or `--n-cells-per-spot-path` followed by the file location.

### Alternative distance metric
By default, CytoSPACE uses Pearson correlation to compare cell and spot transcriptomes. Users can choose to use Spearman correlation or Euclidean distance instead by passing `-dm Spearman_correlation` or `-dm Euclidean` respectively with the function call. 

### Setting a new random seed
While the CytoSPACE algorithm is mostly deterministic, the initial step of sampling cells to be mapped is done at random. To provide an alternative random seed resulting in a different random sampling of cells, users can pass `-se` followed by the desired (integer) seed with the function call. The default random seed for CytoSPACE is 1.

### Alternative handling of sampling
CytoSPACE starts by creating a pool of cells that matches what is expected within the ST data. By default, this is done by resampling single cells to achieve the overall cell type fractions and total cell numbers estimated in the tissue. We recommend that CytoSPACE be run with this default setting for all real data analyses. However, we provide an additional option to generate new "place-holder" cells by sampling from the distribution of gene counts within each cell type instead, and used this option for ensuring uniqueness of mapped cells for benchmarking on simulated data. To run CytoSPACE with this alternative mode, users can pass `-sam place_holders` with the function call. 

### Method extension: mapping quality
While CytoSPACE's formulation as a linear assignment problem guarantees an optimal solution given its cost function, there is no underlying probabilistic framework for estimating mapping uncertainty. One possibility is to determine whether a given cell type belongs to a given spot after mapping - that is, whether a spot contains at least one cell of the same cell type. Notably, this does not distinguish between cells of the same cell type for quality of fit. As such a protocol provides some measure of mapping quality, albeit incomplete, we provide a helper script that implements this via a support vector machine that produces and trains on pseudo-bulks generated from the input scRNA-seq data. This script, `uncertainty_quantification.R`, takes as input the path to the ST dataset count matrix file, the scRNA-seq count matrix file, and the CytoSPACE output file `assigned_locations.csv`, and returns an appended output file with confidence scores in `assigned_locationswConfidenceScores.csv`. The command to run this script following a completed CytoSPACE run is as follows: 
 ```bash
 Rscript uncertainty_quantification.R /path/to/ST_geneexpression /path/to/scRNA_geneexpression /path/to/assigned_locations.csv
```
For interpreting confidence scores, we recommend a cutoff of 0.1, with higher scores indicating increased confidence that a spot contains at least one cell of the same cell type.

### Method extension: single cell ST data
While designed for Visium-type data in which most spots contain RNA from multiple cells, CytoSPACE can also be used with single-cell resolution spatial data such as <a href="https://vizgen.com/resources/meet-the-merscope-platform/" target="_blank">Vizgen's MERSCOPE platform</a>. We expect this extension to be useful for reducing noise and expanding transcriptome coverage of each cell in the ST data. For this single-cell resolution mode, CytoSPACE partitions the ST data into smaller chunks and utilizes multiple CPU cores to assign down-sampled versions of the reference scRNA-seq data to these regions.

To run CytoSPACE with single-cell resolution spatial data:
 ```bash
 cytospace --single-cell 
    --scRNA-path /path/to/scRNA_geneexpression
    --cell-type-path /path/to/scRNA_celllabels
    --st-path /path/to/ST_geneexpression
    --coordinates-path /path/to/ST_coordinates
    --cell-type-fraction-estimation-path path/to/cellfracestimates
    --number-of-processors NUMBER_OF_PROCESSORS
    --number-of-selected-cells  NUMBER_OF_SELECTED_CELLS
    --number-of-selected-spots NUMBER_OF_SELECTED_SPOTS
```
Or with more condensed parameter names: 
 ```bash
 cytospace -sc
    -sp /path/to/scRNA_geneexpression
    -ctp /path/to/scRNA_celllabels
    -stp /path/to/ST_geneexpression
    -cp /path/to/ST_coordinates
    -ctfep path/to/cellfracestimates
    -nop NUMBER_OF_PROCESSORS
    -nosc NUMBER_OF_SELECTED_CELLS
    -noss NUMBER_OF_SELECTED_SPOTS
```
where `NUMBER_OF_PROCESSORS` denotes the number of cores to use, `NUMBER_OF_SELECTED_CELLS` denotes the number of cells in each dowsampled version, and `NUMBER_OF_SELECTED_SPOTS` denotes the size of each ST region. We generally recommend `-nosc 10000 -noss 10000`.

A zip file of example single cell inputs is available to download from Google Drive <a href="https://drive.google.com/file/d/10fhxjCn-VfPPurrI-RE8lbs6NCPqfGXY/view?usp=sharing" target="_blank">here</a>. To download from the command line using `gdown`:
   ```bash
   gdown --fuzzy https://drive.google.com/file/d/10fhxjCn-VfPPurrI-RE8lbs6NCPqfGXY/view?usp=sharing
   unzip single_cell_example_data.zip
   ```

To run CytoSPACE with this example dataset, run the following command from the location of the unzipped inputs and with your CytoSPACE conda environment active:
 ```bash
 cytospace -sc -sp brca_scRNA_GEP.txt -ctp brca_scRNA_celllabels.txt -stp single_cell_ST_example_GEP.txt -cp single_cell_ST_example_coors.txt -ctfep single_cell_ST_example_fractions.txt -nop 1 -nosc 10000 -noss 10000
```

## Extended usage details
```
usage: cytospace [-h] -sp SCRNA_PATH -ctp CELL_TYPE_PATH -stp ST_PATH -cp
                 COORDINATES_PATH -ctfep CELL_TYPE_FRACTION_ESTIMATION_PATH
                 [-ncpsp N_CELLS_PER_SPOT_PATH] [-o OUTPUT_FOLDER]
                 [-op OUTPUT_PREFIX] [-d DELIMITER]
                 [-sm {lapjv,lapjv_compat,lap_CSPR}]
                 [-sam {duplicates,place_holders}] [-mcn MEAN_CELL_NUMBERS]
                 [-se SEED]
                 [-dm {Pearson_correlation,Spearman_correlation,Cosine,Euclidean}]
                 [-nosc NUMBER_OF_SELECTED_CELLS]
                 [-noss NUMBER_OF_SELECTED_SPOTS] [-nop NUMBER_OF_PROCESSORS]
                 [-sc] [-p] [-nr NUM_ROW] [-nc NUM_COLUMN] [-r] [-nv]
                 [-rd ROTATION_DEGREES] [-ss SPOT_SIZE] [-pm PLOT_MARKER]

CytoSPACE is a computational strategy for assigning single-cell transcriptomes
to in situ spatial transcriptomics (ST) data. Our method solves single
cell/spot assignment by minimizing a correlation-based cost function through a
linear programming-based optimization routine.

optional arguments:
  -h, --help            show this help message and exit
  -ncpsp N_CELLS_PER_SPOT_PATH, --n-cells-per-spot-path N_CELLS_PER_SPOT_PATH
                        Path to number of cells per ST spot file
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        Relative path to the output folder
  -op OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Prefix of results stored in the 'output_folder'
  -d DELIMITER, --delimiter DELIMITER
                        Set delimiter of the input files, default '\t'
  -sm {lapjv,lapjv_compat,lap_CSPR}, --solver-method {lapjv,lapjv_compat,lap_CSPR}
                        Which solver to use for the linear assignment problem,
                        default 'lapjv'
  -sam {duplicates,place_holders}, --sampling-method {duplicates,place_holders}
                        Which unerlying method to use for dealing with
                        duplicated cells, default 'duplicates'
  -mcn MEAN_CELL_NUMBERS, --mean-cell-numbers MEAN_CELL_NUMBERS
                        Mean number of cells per spot, default 5 (appropriate
                        for Visium). If analyzing legacy spatial
                        transcriptomics data, set to 20
  -se SEED, --seed SEED
                        Set seed for random generators, default 1
  -dm {Pearson_correlation,Spearman_correlation,Cosine,Euclidean}, --distance-metric {Pearson_correlation,Spearman_correlation,Cosine,Euclidean}
                        Which distance metric to use for the cost matrix,
                        default 'Pearson_correlation'
  -nosc NUMBER_OF_SELECTED_CELLS, --number-of-selected-cells NUMBER_OF_SELECTED_CELLS
                        Number of selected cells from scRNA-seq data used in
                        each iteration
  -noss NUMBER_OF_SELECTED_SPOTS, --number-of-selected-spots NUMBER_OF_SELECTED_SPOTS
                        Number of selected spots from ST data used in each
                        iteration
  -nop NUMBER_OF_PROCESSORS, --number-of-processors NUMBER_OF_PROCESSORS
                        Number of processors used for the analysis
  -sc, --single-cell    Use single-cell spatial approach or not
  -p, --plot-off        Turn create plots on/off
  -nr NUM_ROW, --num-row NUM_ROW
                        Number of rows in pdf figure
  -nc NUM_COLUMN, --num-column NUM_COLUMN
                        Number of columns in pdf figure
  -r, --rotation-flag   Rotate plot
  -nv, --plot-nonvisium
                        Plot with custom slide dimensions
  -rd ROTATION_DEGREES, --rotation-degrees ROTATION_DEGREES
                        Rotation on plot
  -ss SPOT_SIZE, --spot-size SPOT_SIZE
                        Set size of ST spots
  -pm PLOT_MARKER, --plot-marker PLOT_MARKER
                        Shape of ST spots

Required arguments:
  -sp SCRNA_PATH, --scRNA-path SCRNA_PATH
                        Path to scRNA-Seq data
  -ctp CELL_TYPE_PATH, --cell-type-path CELL_TYPE_PATH
                        Path to cell type labels
  -stp ST_PATH, --st-path ST_PATH
                        Path to spatial transcriptomics data (expressions)
  -cp COORDINATES_PATH, --coordinates-path COORDINATES_PATH
                        Path to transcriptomics data (coordinates)
  -ctfep CELL_TYPE_FRACTION_ESTIMATION_PATH, --cell-type-fraction-estimation-path CELL_TYPE_FRACTION_ESTIMATION_PATH
                        Path to cell type fraction file
```

You can see this list of variables and default values for running CytoSPACE from the commmand line as well at any time by calling `cytospace` along with the `-h` or 
`--help` flag, i.e., `cytospace -h`.

## CytoSPACE Solver options
1. `lapjv` __(Recommended for most systems)__    By default, CytoSPACE calls the `lapjv` solver from package `lapjv`. This solver is a fast implementation of the Jonker-Volgenant shortest augmenting path assignment algorithm and returns a globally optimal solution given the objective function as defined in our paper [cite]. As noted above, however, this package is not supported on all systems as it achieves its speedup through use of AVX2 instructions. This solver will be selected by default and can be specified explicitly by passing arguments `--solver-method lapjv` or `-sm lapjv` to `cytospace`.
2. `lap_CSPR` __(Recommended for systems not supporting `lapjv`)__    A second solver option is the `linear_assignment` method from the `ortools` package. This solver uses a different method than the first and third options, an assignment algorithm called the cost scaling push relabel method. This algorithm approximates assignment costs to integer values and loses some numerical precision in doing so. Therefore, while it returns a globally optimal solution __after approximation__ given the objective function defined in the paper, it will return similar but generally not identical results to the first two methods. This solver has a similar running time to the first option and is a good option for systems not supporting the `lapjv` package. This solver can be selected by passing arguments `--solver-method lap_CSPR` or `-sm lap_CSPR` to `cytospace`.
3. `lapjv_compat`   A third solver option implements the `lapjv` solver from package `lap`. Like the first option `lapjv`, this solver also implements the Jonker-Volgenant shortest augmenting path assignment algorithm to return the same globally optimal solution given the objective function defined in the paper. Furthermore, it is broadly supported and should work on all standard operating systems. However, it takes 3-4 times as long to run as the first solver option, the `lapjv` solver from the `lapjv` package, so we only recommend it for systems that do not support the first option. This solver can be selected by passing arguments `--solver-method lapjv_compat` or `-sm lapjv_compat` to `cytospace`.

## Updating local installations
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

## Authors
CytoSPACE was developed in the <a href="https://anlab.stanford.edu/" target="_blank">Newman Lab</a> by

* Milad R. Vahid (miladrv)
* Erin L. Brown (erinlbrown)
* Chloé B. Steen (cbsteen)
* Aaron M. Newman (aaronmnewman)

## Contact
If you have any questions, please contact the CytoSPACE team at cytospaceteam@gmail.com.

## License
CytoSPACE is licensed under the GNU GPL, version 3 or (at your option) any
later version.
CytoSPACE is Copyright (2022-) by the authors.

## Citation
If you use CytoSPACE, please cite:  
  
*Robust alignment of single-cell and spatial transcriptomes with CytoSPACE* (bioxRiv 2022)  
Milad R. Vahid*, Erin L. Brown*,  Chloé B. Steen*,  Minji Kang,  Andrew J. Gentles,  Aaron M. Newman.  
doi: https://doi.org/10.1101/2022.05.20.488356
