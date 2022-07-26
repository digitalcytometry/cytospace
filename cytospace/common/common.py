import numpy as np
import pandas as pd
from pathlib import Path
import warnings
import datatable as dt
from numpy.linalg import norm

def read_file(file_path):
    # Read file
    try:
        file_delim = "," if file_path.endswith(".csv") else "\t"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=pd.errors.ParserWarning)
            #file_data = pd.read_csv(file_path, header=0, index_col=0, delimiter=file_delim)
            file_data = dt.fread(file_path)
            rownames = file_data[:,0].to_pandas().values.flatten()
            file_data = file_data[:,1:file_data.ncols].to_pandas()
            file_data = file_data.set_index(rownames)
   
    except Exception as e:
        raise IOError("Make sure you provided the correct path to input files. "
                      "The following input file formats are supported: .csv with comma ',' as "
                      "delimiter, .txt or .tsv with tab '\\t' as delimiter.")

    return file_data


def normalize_data(data):
    data = np.nan_to_num(data)
    data = 10**6 * (data / np.sum(data, axis=0, dtype=float))
    data_log2 = np.log2(data + 1)
    data_log2 = np.nan_to_num(data_log2)

    return data_log2


def check_paths(output_folder, output_prefix):
    # Create relative path
    output_path = Path().cwd() / output_folder

    # Make sure that the folder exists
    output_path.mkdir(parents=True, exist_ok=True)

    if Path(output_path / f"{output_prefix}assigned_locations.csv").exists():
        print("\033[91mWARNING\033[0m: Running this will overwrite previous results, choose a new"
              " 'output_folder' or 'output_prefix'")

    return output_path


def matrix_correlation_pearson(v1, v2):
    if v1.shape[0] != v2.shape[0]:
        raise ValueError("The two ,atrixes v1 and v2 have to have equal dimensions")

    n = v1.shape[0]
    sums = np.multiply.outer(v2.sum(0), v1.sum(0))
    stds = np.multiply.outer(v2.std(0), v1.std(0))
    correlation = (v2.T.dot(v1) - sums / n) / stds / n

    return correlation


def matrix_correlation_spearman(v1, v2):
    
    if v1.shape[0] != v2.shape[0]:
        raise ValueError("The two ,atrixes v1 and v2 have to have equal dimensions")
        
    v1 = pd.DataFrame(v1).rank().values
    v2 = pd.DataFrame(v2).rank().values
    
    n = v1.shape[0] # v1 and v2 should have the same number of rows
    sums = np.multiply.outer(v2.sum(0), v1.sum(0))
    stds = np.multiply.outer(v2.std(0), v1.std(0))
    correlation = (v2.T.dot(v1) - sums / n) / stds / n
    
    return correlation


def matrix_cosine(v1, v2):
    
    if v1.shape[0] != v2.shape[0]:
        raise ValueError("The two ,atrixes v1 and v2 have to have equal dimensions")
        
    v1 = np.transpose(v1)
    v2 = np.transpose(v2)
       
    dot_product = np.dot(v2,np.transpose(v1)) 
    norm_product = np.outer(norm(v2, axis=1),norm(v1, axis=1))
    cosine = dot_product/norm_product
    
    return cosine
