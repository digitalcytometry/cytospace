<p align="center">
  <img width="300" src="https://github.com/digitalcytometry/cytospace/blob/main/images/CytoSPACE_logo.jpeg">
</p>

<h1> <p align="center">
    Robust and rapid alignment of single-cell and spatial transcriptomes
</p> </h1>

**CytoSPACE** is a novel computational tool for assigning single-cell transcriptomes to in situ spatial transcriptomics (ST) data. Our method solves single cell/spot assignment by minimizing a correlation-based cost function through a linear programming-based optimization routine. 

<p align="center">
  <img src="https://github.com/digitalcytometry/cytospace/blob/main/images/CytoSPACE_overview.png" width="900"> 
</p>

The key innovations of our method are:

- Unlike conventional methods which calculate cell type decompositions by spot, CytoSPACE yields a reconstructed tissue specimen with both high gene coverage and spatially-resolved scRNA-seq data suitable for downstream analysis.
- CytoSPACE is highly robust to noise, and due to its implementation of cell-to-spot assignment via constrained convex optimization, returns globally optimal cell-to-spot assignments. (See the paper for full details.)
- Unlike other methods which generally operate on pre-selected marker genes or on a shared embedding space (the latter of which can erase true biological variation), CytoSPACE uses the full transcriptome without the need for batch correction, helping it retain sensitivity to subtle cell states.
- CytoSPACE is quick and simple to execute. It runs in minutes even with a single CPU on a personal laptop and requires no hyperparameter tuning or gene/feature selection.

## Installation instructions
1. Install <a href="https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html" target="_blank">Miniconda</a> if not already available.

2. Clone this repository:
```bash
git clone https://github.com/digitalcytometry/cytospace
```

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

7. (Recommended) Install package `lapjv` by executing:
```bash
pip install lapjv
```
We highly recommend you install this package, which provides a fast implementation of the default core optimization algorithm within CytoSPACE. However, some systems may not accommodate it as it requires CPU support for AVX2 instructions. To determine if your system supports this package, it is generally easiest to simply attempt to install it as above. If it installs without problems, your system will support it! If you run into an error, it is likely your system does not support it, and you can simply use one of the other options we have provided. See __Solver options__ below for details. 

## File format
CytoSPACE requires 5 files as input. All files should be provided in tab or comma-delimited tabular input format (saved as .txt or .csv, respectively) with no double quotations. Further formatting details for each input file are specified below:

1. __A scRNA-seq gene expression file:__
- The matrix must be genes (rows) by cells (columns).
- The first row must contain the single cell IDs and the first column must contain the gene names.
- The first column (gene names) must have a header.
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
- The first column (gene names) must have a header.
- The gene expression data should be represented as non-normalized counts. 
<img src="https://github.com/digitalcytometry/cytospace/blob/main/images/STdatafile.png" width="800"> 

4. __A spatial transcriptomics coordinates file:__
- A table consisting of 3 columns, where the first column contains the ST spot IDs corresponding to and in the same order as the columns of the ST gene expression matrix, and column 2 and 3 contain the row and column indices of the spatial transcriptomics data, respectively. 
- The columns must have a header. 
<img src="https://github.com/digitalcytometry/cytospace/blob/main/images/STcoordfile.png" width="250"> 

5. __A file with cell type fraction estimates, obtained from the `R` script `get_cellfracs_seuratv3.R`.__ 
- A table consisting of 2 rows, where the first row is the cell type labels, and the second row is the cell fractions of each cell type represented as proportions between 0 and 1. The first column is the row names. 
- For further details on running `get_cellfracs_seuratv3.R`, see section "__Preprocessing__" below.
<img src="https://github.com/digitalcytometry/cytospace/blob/main/images/cell_type_fractions_file.png" width="800"> 

                                                                                                                 
## File preparation
If you are starting with outputs from Cell Ranger (scRNA-seq from 10x) or Space Ranger (ST from 10x), you can use the `R` script `???` to produce files formatted for CytoSPACE input.

```bash
ADD DETAILS
```

Similarly, if you are starting with Seurat objects derived from any source, you can use the `R` function `???` to produce files formatted for CytoSPACE input.

```bash
ADD DETAILS
```

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
2. We use `Seurat v3` for estimating cell fractions, and this installed as part of the CytoSPACE environment. If you want to run other analyses using more recent versions of Seurat after running CytoSPACE, for example Seurvat v4, make sure to first __deactivate the CytoSPACE environment__ once you are done running CytoSPACE. This is done using the command `deactivate cytospace`.

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

For full usage details with additional options, see the __Extended usage details__ section below. 

### Choosing a solver
CytoSPACE provides three solver options. In short, we recommend using the default option `lapjv` if your system supports AVX2 (i.e., if you were able to successfully install it with `pip install lapjv`) and `lap_CSPR` otherwise. No options are required to use the default solver `lapjv`. To use `lap_CSPR` instead, pass the argument `-sm lap_CSPR` to your `cytospace` call. For full solver details, see the __Solver options__ section below.


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

## CytoSPACE outputs
CytoSPACE will produce five output files by default.
1. ```FILE NAME``` Heatmaps of cell type assignments within the ST sample. Along with a plot showing the total number of cells mapped to each spot, these will show the number of cells per cell type mapped to each spot.
2. ```FILE NAME``` This file will provide the assigned locations of each single cell mapped to ST spots. As some cells may be mapped to multiple locations depending on the size of the input scRNA-seq set, new cell IDs (`UniqueCID`) are assigned to each cell and given in the first column. The second column includes original cell IDs (`OriginalCID`); the third includes assigned spot IDs (`SpotID`); the fourth and fifth columns respectively include  `row` and `column` indices of the corresponding spots; and then optionally, the sixth and seventh columns include `coord_x` and `coord_y` of the corresponding spots if these details were provided in the initial Coordinates file.
3. ```FILE NAME``` This file gives the raw number of cells of each cell type per spot by `SpotID` as well as the total number of cells assigned to that spot.
4. ```FILE NAME``` This file gives the fractional abundance of cell types assigned to each spot by `SpotID`.
5. ```log.txt``` This file contains a log of CytoSPACE run parameters and running time.

## Example datasets for running CytoSPACE
For users to test CytoSPACE, we have included files for two example runs:
1. A HER2+ breast cancer scRNA-seq atlas by Wu et al. (<a href="https://www.nature.com/articles/s41588-021-00911-1" target="_blank">Nature Genetics, 2021</a>) and a HER2+ breast cancer FFPE specimen profiled by the Visium platform (<a href="https://www.10xgenomics.com/resources/datasets/human-breast-cancer-ductal-carcinoma-in-situ-invasive-carcinoma-ffpe-1-standard-1-3-0" target="_blank">10x Genomics</a>). Default parameters were selected with Visium samples in mind and are appropriate here.
2. A melanoma scRNA-seq atlas by Tirosh et al. (<a href="https://www.science.org/doi/10.1126/science.aad0501?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed" target="_blank">Science, 2016</a>), and a melanoma specimen profiled by the legacy ST platform (Thrane et al, <a href="https://aacrjournals.org/cancerres/article/78/20/5970/631815/Spatially-Resolved-Transcriptomics-Enables" target="_blank">Cancer Research, 2018</a>). As the legacy ST platform has larger spot sizes, we recommend mapping an average of 20 cells per spot, i.e. passing the argument `-mcn 20` to your `cytospace` call.

### Download example datasets
Zip files containing the example datasets can be downloaded from the following links:
1. <a href="https://drive.google.com/file/d/1vAqszYk3-B2vgwkSFMprsUcRBFr-lS2f/view?usp=sharing" target="_blank">Breast cancer</a>
2. <a href="https://drive.google.com/file/d/1oHe4UP2K0kQS9gNFVvZtjqtJyJK_VeBp/view?usp=sharing" target="_blank">Melanoma</a> 

To download from the command line using `gdown`:
1. Breast cancer
   ```bash
   gdown --fuzzy https://drive.google.com/file/d/1vAqszYk3-B2vgwkSFMprsUcRBFr-lS2f/view?usp=sharing
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
cytospace -sp brca_scRNA_GEP.txt -ctp brca_scRNA_celllabels.txt -stp brca_STdata_GEP.txt -cp brca_STdata_coordinates.txt -ctfep brca_cell_fraction_estimates.txt
```

```bash
cytospace -sp melanoma_scRNA_GEP.txt -ctp melanoma_scRNA_celllabels.txt -stp melanoma_STdata_slide1_GEP.txt -cp melanoma_STdata_slide1_coordinates.txt -ctfep melanoma_cell_fraction_estimates.txt -mcn 20
```

## CytoSPACE output files for example breast cancer data
The main output from a CytoSPACE run is the file named `assigned_locations.csv`, which provides the ST spots to which the single cells have been assigned. 

```
include image here
```

The CytoSPACE results are visualized in heatmaps saved as `plot_cell_type_locations.pdf` showing the distribution of single cells across ST spots for each cell type. Below are the heatmaps produced for the example BRCA data.

<p align="center">
  <img width="800" src="https://github.com/digitalcytometry/cytospace/blob/main/images/BRCA_plot_cell_type_locations.png">
</p>

For comparison, consider the pathologist annotations of this ST sample as provided by 10x:

<p align="center">
  <img width="800" src="https://github.com/digitalcytometry/cytospace/blob/main/images/Visium_FFPE_Human_Breast_Cancer_Pathologist_Annotations.png">
</p>

The number of cells per spot by cell type as well as in total are provided in the file `cell_type_assignments_by_spot.csv`.
```
include image here
```
Fractional abundances of each cell type are returned in the file `fractional_abundances_by_spot.csv`
```
include image here
```

## Extended usage details
```
usage: cytospace [-h] -sp SCRNA_PATH -ctp CELL_TYPE_PATH -stp ST_PATH -cp COORDINATES_PATH -ctfep
                 CELL_TYPE_FRACTION_ESTIMATION_PATH [-o OUTPUT_FOLDER] [-op OUTPUT_PREFIX] [-d DELIMITER]
                 [-m {shortest_augmenting_path,cost_scaling_push_relabel}] [-sm {lap,lapjv}] [-mcn MEAN_CELL_NUMBERS]
                 [-se SEED] [-p] [-nr NUM_ROW] [-nc NUM_COLUMN] [-r] [-rd ROTATION_DEGREES] [-ss SPOT_SIZE]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        Relative path to the output folder, default 'cytospace_results'
  -op OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Prefix of results stored in the 'output_folder', default ''
  -d DELIMITER, --delimiter DELIMITER
                        Set delimiter of the input files, default ',' (set to '\t' for tab-delimited files)
  -m {shortest_augmenting_path,cost_scaling_push_relabel}, --method {shortest_augmenting_path,cost_scaling_push_relabel}
                        Method for computing the linear assignment sum, default 'shortest_augmenting_path'
  -sm {lap,lapjv}, --solver-method {lap,lapjv}
                        Which solver to use for the linear assignment problem when setting
                        'method'='shortest_augmenting_path', default 'lapjv'
  -mcn MEAN_CELL_NUMBERS, --mean-cell-numbers MEAN_CELL_NUMBERS
                        Mean number of cells per spot, default 5 (appropriate for most Visium -- if analyzing legacy spatial
                        transcriptomics data, other options e.g. 20 may be more appropriate)
  -se SEED, --seed SEED
                        Set seed for random generators, default 1
  -p, --plot-off        Turn create plots on/off, default False (set to True if you do not want plots)
  -nr NUM_ROW, --num-row NUM_ROW
                        Number of rows in pdf figure, default 4
  -nc NUM_COLUMN, --num-column NUM_COLUMN
                        Number of coulmns in pdf figure, default 4
  -r, --rotation-flag   Rotate plot, default True (appropriate for Visium data with row and column indices provided)
  -rd ROTATION_DEGREES, --rotation-degrees ROTATION_DEGREES
                        Rotation on plot, default 270 (appropriate for Visium data with row and column indices provided)
  -ss SPOT_SIZE, --spot-size SPOT_SIZE
                        Set size of ST spots, default 155 (appropriate for standard Visium data)

Required arguments:
  -sp SCRNA_PATH, --scRNA-path SCRNA_PATH
                        Path to scRNA-Seq data, which should be a 
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

## Solver options
1. `lapjv` __(Recommended for most systems)__    By default, CytoSPACE calls the `lapjv` solver from package `lapjv`. This solver is a fast implementation of the Jonker-Volgenant shortest augmenting path assignment algorithm and returns a globally optimal solution given the objective function as defined in our paper [cite]. As noted above, however, this package is not supported on all systems as it achieves its speedup through use of AVX2 instructions. This solver will be selected by default and can be specified explicitly by passing arguments `--solver-method lapjv` or `-sm lapjv` to `cytospace`.
2. `lap_CSPR` __(Recommended for systems not supporting `lapjv`)__    A second solver option is the `linear_assignment` method from the `ortools` package. This solver uses a different method than the first and third options, an assignment algorithm called the cost scaling push relabel method. This algorithm approximates assignment costs to integer values and loses some numerical precision in doing so. Therefore, while it returns a globally optimal solution __after approximation__ given the objective function defined in the paper, it will return similar but generally not identical results to the first two methods. This solver has a similar running time to the first option and is a good option for systems not supporting the `lapjv` package. This solver can be selected by passing arguments `--solver-method lap_CSPR` or `-sm lap_CSPR` to `cytospace`.
3. `lapjv_compat`   A third solver option implements the `lapjv` solver from package `lap`. Like the first option `lapjv`, this solver also implements the Jonker-Volgenant shortest augmenting path assignment algorithm to return the same globally optimal solution given the objective function defined in the paper. Furthermore, it is broadly supported and should work on all standard operating systems. However, it takes 3-4 times as long to run as the first solver option, the `lapjv` solver from the `lapjv` package, so we only recommend it for systems that do not support the first option. This solver can be selected by passing arguments `--solver-method lapjv_compat` or `-sm lapjv_compat` to `cytospace`.

## Updating local installations
To update your local installation of CytoSPACE following updates of this GitHub repository, navigate to your `cytospace` directory and execute the following commands:
```bash
git pull
conda env update --name cytospace --file environment_withoutlapjv.yml
pip install .
```
If you have made local updates to your version of the CytoSPACE source code, you should execute 
```bash
pip install .
``` 
once more before running. 

## Authors
CytoSPACE was developed by

* Milad R. Vahid (miladrv)
* Erin L. Brown (erinlbrown)
* Chloé B. Steen (cbsteen)
* Aaron M. Newman (aaronmnewman)

## License
CytoSPACE is licensed under the GNU GPL, version 3 or (at your option) any
later version.
CytoSPACE is Copyright (2022-) by the authors.

## Citation
If you use CytoSPACE, please cite:

    Coming soon
