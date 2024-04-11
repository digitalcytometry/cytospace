import numpy as np
import pandas as pd
import scanpy as sc

import datatable as dt
import scipy.io
import tarfile

import os
import subprocess
from pathlib import Path

import warnings
import pkg_resources

def read_file(file_path):
    # Read file
    is_sparse = (file_path.endswith(".mtx") or file_path.endswith(".mtx.gz"))
    if is_sparse:
        try:
            if not os.path.isfile(file_path):
                raise IOError("Cannot locate file: {}".format(file_path))

            # Lists of possible names and extensions depending on CellRanger version
            possible_gene_names = ["genes", "features"]
            possible_cell_names = ["cells", "barcodes"]
            
            possible_extensions = [".tsv", ".csv", ".tsv.gz", ".csv.gz"]

            # Function to find and read file based on possible names and extensions
            def find_file(base_path, possible_names, possible_extensions):
                for name in possible_names:
                    for extension in possible_extensions:
                        full_path = f"{base_path}{name}{extension}"
                        if os.path.isfile(full_path):
                            return full_path, extension
                raise IOError(f"Required files not found for base path: {base_path}")

            base_file_path = os.path.dirname(file_path) + os.path.sep

            genes_path, genes_extension = find_file(base_file_path, possible_gene_names, possible_extensions)
            cells_path, cells_extension = find_file(base_file_path, possible_cell_names, possible_extensions)
            genes = pd.read_csv(genes_path, sep='\t' if '.tsv' in genes_extension else ',', header=None).iloc[:, 0].to_numpy()
            cells = pd.read_csv(cells_path, sep='\t' if '.tsv' in cells_extension else ',', header=None).iloc[:, 0].to_numpy()

            sparse_mtx = scipy.io.mmread(file_path)

            if not ((sparse_mtx.shape[0] == len(genes)) and (sparse_mtx.shape[1] == len(cells))):
                raise IOError("""
                    The dimensions of the provided sparse matrix does not match the corresponding gene and cell lists.
                    The files for gene and cell lists must have one entry on each line, without a header line.
                    Please check the following files: {}, {}.
                """.format(genes_path, cells_path))
            
            file_data = pd.DataFrame.sparse.from_spmatrix(sparse_mtx, index=genes, columns=cells).sparse.to_dense()

        except Exception as e:
            print("Error encountered while reading in sparse matrix as input: {}".format(file_path))
            raise
            
    else:
        try:
            file_delim = "," if file_path.lower().endswith(".csv") else "\t"

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=pd.errors.ParserWarning)
                file_data = dt.fread(file_path, header=True)
                colnames = pd.read_csv(file_path, sep=file_delim, nrows=1, index_col=0).columns
                rownames = file_data[:, 0].to_pandas().values.flatten()
                file_data = file_data[:, 1:].to_pandas()
                file_data.index = rownames
                file_data.columns = colnames

        except Exception as e:
            print("Error encountered while reading input: {}".format(file_path))
            print("Please make sure that you provided the correct path to the input files.",
                  "The following input file formats are supported:",
                  ".csv with comma ',' as delimiter,",
                  "and .txt or .tsv with tab '\\t' as delimiter.")
            raise

    return file_data


def read_visium(file_name,out_dir):
    sr_dir = os.path.join(out_dir, 'SpaceRanger')

    file = tarfile.open(file_name,'r')
    file.extractall(sr_dir)
    file.close()

    count_dir = os.listdir(sr_dir)
    count_dir = [d for d in count_dir if os.path.isdir(os.path.join(sr_dir, d))][0]
    count_dir = os.path.join(sr_dir, count_dir)

    fin_count = os.listdir(count_dir)
    fin_count = [fin for fin in fin_count if fin.endswith('.h5')]
    fin_count = [fin for fin in fin_count if not fin.startswith('.')][0]

    # check for the presence of "tissue_positions_list.csv" file
    # if not present, read "tissue_positions.csv" file and write as "tissue_positions_list.csv"
    if 'tissue_positions_list.csv' not in os.listdir(count_dir+"/spatial"):
        if 'tissue_positions.csv' in os.listdir(count_dir+"/spatial"):
            df = pd.read_csv(count_dir+"/spatial/tissue_positions.csv")
            df.to_csv(count_dir+"/spatial/tissue_positions_list.csv",index=False)
    
    visium_adata = sc.read_visium(count_dir,count_file=fin_count)
    visium_adata.var_names_make_unique(join='.')

    genes = visium_adata.var_names
    spots = visium_adata.obs_names
    df_expression = pd.DataFrame.sparse.from_spmatrix(visium_adata.X,index=spots,columns=genes).transpose().sparse.to_dense()

    df_expression.index.name = 'GENES'

    df_coords = visium_adata.obs[['array_row', 'array_col']]
    df_coords.columns = ['row','col']
    df_coords = df_coords.loc[df_expression.columns,:]
    df_coords.index.name = 'SpotID'

    return df_expression, df_coords

def estimate_cell_type_fractions(scRNA_path, cell_type_path, st_path, output_folder, output_prefix):
    # Get path to R script
    run_script = pkg_resources.resource_filename("cytospace", "get_cellfracs_seuratv3.R")

    # Windows paths
    scRNA_path = scRNA_path.replace("\\", "\\\\")
    cell_type_path = cell_type_path.replace("\\", "\\\\")
    st_path = st_path.replace("\\", "\\\\")
    output_folder = output_folder.replace("\\", "\\\\")

    # Run command
    run_args = ["Rscript", run_script,
                "--scrna-path", scRNA_path, "--ct-path", cell_type_path, "--st-path", st_path,
                "--outdir", output_folder, "--prefix", output_prefix]
    out = subprocess.run(run_args, check=True)

    frac_file = os.path.join(output_folder, ''.join([output_prefix, 'Seurat_weights.txt']))
    return frac_file

def normalize_data(data):
    data = np.nan_to_num(data).astype(float)
    data *= 10**6 / np.sum(data, axis=0, dtype=float)
    np.log2(data + 1, out=data)
    np.nan_to_num(data, copy=False)
    return data

def downsample(data_df, target_count):
    """
    Downsamples dataset to target_count transcripts per cell.

    Parameters :
        data_df (pd.DataFrame) :
            data to downsample, with rows as genes and columns as cells.

    Returns :
        downsampled_df (pd.DataFrame) :
            downsampled data, where the sum of each column is target_count
            (or lower, for columns whose sum was originally lower than target_count).
    """
    def downsample_cell(sr, target_tr_count):
        if sr.sum() <= target_tr_count:
            return sr

        genes, counts = np.unique(np.random.choice(np.repeat(sr.index, sr.to_numpy()), target_tr_count), return_counts=True)
        downsampled = pd.Series(counts, index=genes).reindex(sr.index, fill_value=0)

        return downsampled
    
    downsampled_df = data_df.apply(lambda k: downsample_cell(k, target_count), axis=0)

    return downsampled_df


def check_paths(output_folder, output_prefix):
    # Create relative path
    output_path = os.path.join(os.getcwd(), output_folder)

    # Make sure that the folder exists
    Path(output_path).mkdir(parents=True, exist_ok=True)

    if os.path.exists(os.path.join(output_path, f"{output_prefix}assigned_locations.csv")):
        print("\033[91mWARNING\033[0m: Running this will overwrite previous results, choose a new"
              " 'output_folder' or 'output_prefix'")

    return output_path


def matrix_correlation_pearson(v1, v2):
    if v1.shape[0] != v2.shape[0]:
        raise ValueError("The two matrices v1 and v2 must have equal dimensions; ST and scRNA data must have the same genes")

    n = v1.shape[0]
    sums = np.multiply.outer(v2.sum(0), v1.sum(0))
    stds = np.multiply.outer(v2.std(0), v1.std(0))
    correlation = (v2.T.dot(v1) - sums / n) / stds / n

    return correlation


def matrix_correlation_spearman(v1, v2):
    
    if v1.shape[0] != v2.shape[0]:
        raise ValueError("The two matrices v1 and v2 must have equal dimensions; ST and scRNA data must have the same genes")
        
    v1 = pd.DataFrame(v1).rank().values
    v2 = pd.DataFrame(v2).rank().values
    
    n = v1.shape[0] # v1 and v2 should have the same number of rows
    sums = np.multiply.outer(v2.sum(0), v1.sum(0))
    stds = np.multiply.outer(v2.std(0), v1.std(0))
    correlation = (v2.T.dot(v1) - sums / n) / stds / n
    
    return correlation
