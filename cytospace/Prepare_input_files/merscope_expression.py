"""
merscope_expression.py

Defines the MerscopeExpression class, to be used in process_merscope.py.
"""

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import h5py
import matplotlib.path as mplpath

import os
import sys
from pathlib import Path

import time


# Suffix for CytoSPACE input files
CYTOSPACE_ST_EXPRESSION_SUFFIX = 'ST_expressions_cytospace.tsv'
CYTOSPACE_ST_CELLTYPE_SUFFIX = 'ST_annotations_cytospace.tsv'
CYTOSPACE_ST_COORDINATE_SUFFIX = 'ST_coordinates_cytospace.tsv'
CYTOSPACE_SC_EXPRESSION_SUFFIX = 'scRNA_expressions_cytospace.tsv'
CYTOSPACE_SC_CELLTYPE_SUFFIX = 'scRNA_annotations_cytospace.tsv'


class MerscopeExpression(object):
    def __init__(self, sample, data_dir, save_dir=None, save_intermediates=False):
        self.sample = sample        # str, e.g., "HumanBreastCancerPatient1"
        self.data_dir = data_dir    # str, e.g., "/path/to/raw/data/HumanBreastCancerPatient1/"

        self.save_dir = None        # str, e.g., "/path/to/save/outputs/"
        if save_dir is None:
            print('Setting save_dir to the current directory...', flush=True)
            self.save_dir = os.getcwd()
        else:
            Path(save_dir).mkdir(parents=True, exist_ok=True)
            self.save_dir = save_dir
        self.save_intermediates = save_intermediates # bool

        self.raw_expr = None        # pd.DataFrame, cell x gene matrix output
        self.gene_list = None       # List[str]
        self.metadata = None        # pd.DataFrame, with cell ID as index and columns "X", "Y"

        self.annotations = None     # pd.DataFrame, with cell ID as index and column name "cell_type"

        self.z_exprs = {}           # Dict{z_idx(int) : cell x gene matrix(pd.DataFrame)}

        self.z_peri_idx = None      # int
        self.z_mid_idx = None       # int
        self.cell_types = None      # List[str]
        self.sanitized_expr = None  # pd.DataFrame, cell x gene matrix

    def read_raw_data(self, raw_expr_file="cell_by_gene.csv", sep=',', is_abs_path=False, remove_blank_genes=True):
        """
        Reads in the raw count matrix (self.raw_expr).

        raw_expr_file (str) :
            Name of CSV file to read count matrix (cell x gene) from.
            By default, reads cell_by_gene.csv from self.data_dir.
        sep         (str) : Separator character; ',' by default.
        is_abs_path (bool): raw_expr_file is treated as absolute path if True,
                            or as relative path to self.data_dir if False.
                            False by default.
        remove_blank_genes (bool) : Removes "Blank-nn" genes from expression matrix if True.
                                    True by default.
        
        Returns None; sets self.raw_expr, self.gene_list.
        """
        raw_expr_path = raw_expr_file if is_abs_path else os.path.join(self.data_dir, raw_expr_file) 

        self.raw_expr = pd.read_csv(raw_expr_path, index_col=0, sep=sep, dtype=int)
        self.gene_list = self.raw_expr.columns.tolist()

        if remove_blank_genes:
            # remove blank genes
            self.gene_list = [gene for gene in self.gene_list if not gene.startswith('Blank')]
            self.raw_expr = self.raw_expr.loc[:, self.gene_list]

    def read_metadata(self, metadata_file="cell_metadata.csv", col_X='center_x', col_Y='center_y', sep=',', is_abs_path=False):
        """
        Reads in the metadata (self.metadata).

        metadata_file (str) :
            Name of CSV file to read metadata from.
            By default, reads cell_metadata.csv from self.data_dir.
        col_X, col_Y (str) : Column names representing xy-coordinates for each cell.
        sep          (str) : Separator character; ',' by default.
        is_abs_path (bool) : metadata_file is treated as absolute path if True,
                             or as relative path to self.data_dir if False.
                             False by default.
        
        Returns None; sets self.metadata.
        """
        metadata_path = metadata_file if is_abs_path else os.path.join(self.data_dir, metadata_file) 

        self.metadata = pd.read_csv(metadata_path, sep=sep, index_col=0)
        self.metadata.rename(columns={col_X : 'X', col_Y : 'Y'}, inplace=True)

    def set_annotations(self, annotation_file, col_name, sep=',', is_abs_path=False):
        """
        Initializes annotations from a file.

        annotation_file (str) :
            Name of CSV file to read annotations from, with the row names being cell IDs.
        col_name    (str) : Name of column listing the cell types.
        sep         (str) : Separator character; ',' by default.
        is_abs_path (bool): annotation_file is treated as absolute path if True,
                            or as relative path to self.data_dir if False.
                            False by default.

        Returns None; sets self.annotations.
        """
        annotation_path = annotation_file if is_abs_path else os.path.join(self.data_dir, annotation_file)

        data = pd.read_csv(annotation_path, index_col = 0, sep = sep)
        annotations = data[[col_name]]
        annotations.columns = ['cell_type']

        self.annotations = annotations
        if self.sanitized_expr is None:
            self.cell_types = annotations.cell_type.unique().tolist()

    def _get_tr_count(self, transcripts_fov, z_idx, cell_boundaries_abs_path, verbose):
        """
        For a transcript subset inside a single field-of-view, aggregates the transcript count for each gene on the specified z-index.
        Designed to be called inside self.get_zplanes().
        """
        t0 = time.time()
        # identify relevant FOV: transcripts dataframe should have been grouped by FOV when calling this function
        if len(transcripts_fov.fov.unique()) > 1:
            raise ValueError('_get_tr_count() error: nonunique values in fov ({})'.format(transcripts_fov.fov.unique()))
        fov = transcripts_fov.fov.unique()[0]
        # read corresponding HDF5 file
        hdf5_abs_path = os.path.join(cell_boundaries_abs_path, 'feature_data_{}.hdf5'.format(fov))
        if not os.path.isfile(hdf5_abs_path):
            print('Warning: file feature_data_{}.hdf5 is missing in {}; skipping...'.format(fov, cell_boundaries_abs_path))
            return pd.DataFrame(columns=self.gene_list, dtype=int)
        # initialize variables
        coord_obj = h5py.File(hdf5_abs_path)
        expr = pd.DataFrame(np.zeros((len(coord_obj['featuredata'].keys()), len(self.gene_list))),
                            index=coord_obj['featuredata'].keys(), columns=self.gene_list, dtype=int) # cell x gene dataframe, initialized with zeros
        # count transcript for each cell
        for cell_id in coord_obj['featuredata'].keys():
            coord_arr = np.array(coord_obj['featuredata'][cell_id]['zIndex_' + str(z_idx)]['p_0']['coordinates'])
            path_coords = mplpath.Path(coord_arr[0], closed=True)
            min_xy = np.min(coord_arr[0], axis=0)
            max_xy = np.max(coord_arr[0], axis=0)
            transcripts_nearby = transcripts_fov.loc[(transcripts_fov['global_x'] >= min_xy[0]) &
                                                        (transcripts_fov['global_x'] <= max_xy[0]) &
                                                        (transcripts_fov['global_y'] >= min_xy[1]) &
                                                        (transcripts_fov['global_y'] <= max_xy[1])]
            transcript_is_in_cell = path_coords.contains_points(transcripts_nearby[['global_x', 'global_y']].values)
            count_by_gene_zplane = transcripts_nearby.loc[transcript_is_in_cell].gene.value_counts().reindex(self.gene_list)
            expr.loc[cell_id] = count_by_gene_zplane
        if verbose:
            print('Completed _get_tr_count() for fov={}: time elapsed {}'.format(fov, time.time() - t0))
        return expr

    def get_zplanes(self, z_idxs,
                    cell_boundaries_dir='cell_boundaries/', transcript_file='detected_transcripts.csv',
                    sep=',', is_abs_path=False, regenerate=False, verbose=False, save_intermediates=None):
        """ 
        Generates expression matrix for a single z-plane.
        Note: This function may take some time (~1 hour) to complete.

        z_idxs (int or List[int]) : Index(es) of z-plane to generate the expression matrix for.
        cell_boundaries_dir (str) : Directory containing HDF5 files for cell polygon coordinates.
                                    By default, reads from "cell_boundaries" under self.data_dir.
        transcript_file     (str) : Name of CSV file listing transcript data.
                                    Requires columns "global_x", "global_y", "global_z", "gene", "fov".
                                    By default, reads from "detected_transcripts.csv" under self.data_dir.
        sep         (str) : Separator character for transcript_file; ',' by default.
        is_abs_path (bool): cell_boundaries_dir and transcript_file are treated as absolute paths if True,
                            or as relative paths to self.data_dir if False.
                            False by default.
        regenerate  (bool): If True, regenerates the matrix even if one already exists in self.z_expr.
        verbose     (bool): Prints progress if True. False by default.
        save_intermediates (bool) :
            Whether to save z-plane count matrix as a CSV file.
            By default, follows self.save_intermediates.
        
        Returns None; adds an entry to self.z_exprs.
        """
        if isinstance(z_idxs, int):
            z_idxs = [z_idxs]
        
        transcripts_full = None
        for z_idx in z_idxs:
            if z_idx in self.z_exprs.keys() and not regenerate:
                print(f"Gene expression matrix for z={z_idx} already exists; rerun get_zplanes() with regenerate=True if you would like to overwrite")
                continue
            
            transcript_path = transcript_file if is_abs_path else os.path.join(self.data_dir, transcript_file)
            cell_boundaries_path = cell_boundaries_dir if is_abs_path else os.path.join(self.data_dir, cell_boundaries_dir)

            if transcripts_full is None:
                transcripts_full = pd.read_csv(transcript_path, sep=sep,
                                                usecols=['global_x', 'global_y', 'global_z', 'gene', 'fov'],
                                                dtype={'global_x':'float32', 'global_y':'float32',
                                                        'global_z':'uint8', 'gene':'category', 'fov':'uint16'})
                transcripts_full = transcripts_full[(transcripts_full['global_z'] == z_idx) &
                                                    (transcripts_full['gene'].isin(self.gene_list))]
            
            if verbose:
                print('Initiating _get_tr_count() for z={} over {} FOVs...'.format(z_idx, len(transcripts_full.fov.unique())))

            res = transcripts_full.groupby('fov').apply(lambda x: self._get_tr_count(x, z_idx, cell_boundaries_path, verbose))
            res.reset_index(level=[0], inplace=True)
            res.drop('fov', axis=1, inplace=True)

            self.z_exprs[z_idx] = res

            if save_intermediates or (save_intermediates is None and self.save_intermediates):
                print('Saving gene expression matrix for z={}...'.format(z_idx))
                self.z_exprs[z_idx].to_csv(os.path.join(self.save_dir, '{}_cell_by_gene_z{}.csv'.format(self.sample, z_idx)))

    def sanitize_matrix(self, z_peri_idx=0, z_mid_idx=3, exclude_celltypes=None, pval_cutoff=0.05, save_intermediates=None):
        """
        Removes expression of potentially contaminating genes for each cell type, identified from distribution across z-planes.

        z_peri_idx    (int) : Index of z-plane to use as the "peripheral" plane. 0 by default.
        z_mid_idx     (int) : Index of z-plane to use as the "center" plane. 3 by default (out of 0-6).
        z_orig_expr_file (List[str]) :
            List of annotated cell types to ignore (thus not removing contamination), if any. None by default.
        pval_cutoff (float) : P-value threshold for differential expression between peripheral and center z-planes. 0.05 by default.
        save_intermediates (bool) :
            Whether to save sanitized counts as a CSV file.
            By default, follows self.save_intermediates.
        
        Returns None; sets self.z_peri_idx, self.z_mid_idx, self.sanitized_expr.
        """
        if self.raw_expr is None:
            raise ValueError('Raw expression matrix has not been correctly added; please first call read_raw_data()')
        if self.annotations is None:
            raise ValueError('Annotations have not been correctly added; please first call set_annotations()')
        if z_peri_idx not in self.z_exprs.keys():
            raise ValueError(f'Gene expression matrix for z={z_peri_idx} has not been generated; please first call get_zplanes() for z={z_peri_idx}')
        if z_mid_idx not in self.z_exprs.keys():
            raise ValueError(f'Gene expression matrix for z={z_mid_idx} has not been generated; please first call get_zplanes() for z={z_mid_idx}')

        self.z_peri_idx = z_peri_idx
        self.z_mid_idx = z_mid_idx
        
        # define peripheral z-plane
        z_edge = self.z_exprs[z_peri_idx]
        # define center z-plane
        z_mid = self.z_exprs[z_mid_idx]

        # initialize list of cell types, excluding some if specified in exclude_celltypes
        celltype_list = np.sort(self.annotations.cell_type.unique())
        celltype_list = np.delete(celltype_list, np.where(np.isin(celltype_list, exclude_celltypes)))
        self.cell_types = celltype_list.tolist()

        # initialize sanitized matrix as the raw matrix
        # note: sanitization is based on cell type, so sanitized matrix will be subsetted to those in celltype_list
        cell_idxs = self.annotations.loc[self.annotations.cell_type.isin(celltype_list)].index
        self.sanitized_expr = self.raw_expr.loc[cell_idxs, :]

        pvals_list = [] # List[List[np.float]], one entry per celltype
        markers = {}    # { celltype(str) : List[gene(str)] }

        # identify potential source genes of contamination
        for celltype in celltype_list:
            pvals_ctype = [] # List[np.float], one entry per gene
            print(celltype, flush=True)
            
            # subset to cells of this cell type
            cells = self.annotations[self.annotations.cell_type == celltype].index

            for gene in self.gene_list:
                # get the expressions of the current gene from the center and the peripheral z-plane
                z_edge_expr = z_edge.loc[cells, gene]
                z_mid_expr = z_mid.loc[cells, gene]
                
                try:
                    # test significance
                    w = wilcoxon(z_mid_expr, z_edge_expr)
                    pval = w.pvalue
                    if pval == 0:
                        pval = sys.float_info.min
                        # if pval == 0.0, it may be smaller than the minimum value (~2e-308) representable by float
                        # assign the minimum value in this case
                    # compute pvalue score (to be plotted later if necessary)
                    pval = -np.log10(pval) * np.sign(z_mid_expr.sum() - z_edge_expr.sum())
                    pvals_ctype.append(pval)
                except Exception as e:
                    # may give error if all counts are zero for either of the z-planes
                    print('Zero warning {}/{}: '.format(celltype, gene) + str(e), flush=True)
                    pvals_ctype.append(np.nan)
                    continue
            
            pvals_ctype = np.array(pvals_ctype)
            markers_ctype = self.gene_list[pvals_ctype > -np.log10(pval_cutoff) * (1)].tolist()
            markers[celltype] = markers_ctype
            pvals_list.append(pvals_ctype)

        # pvals: pd.DataFrame of -log10(p); positive if the expression is higher in z_mid than z_edge, negative if not
        pvals = pd.DataFrame(pvals_list)
        pvals.index = celltype_list
        pvals.columns = self.gene_list

        # aggregate all contamination candidate genes
        markers_all = []
        for celltype, ct_markers in markers.items():
            markers_all.extend(ct_markers)
        markers_all = set(markers_all)

        # identify contamination using a similar method
        contamination_genes = {} # { celltype(str) : List[gene(str)] }
        for celltype in celltype_list:
            # exclude candidates that were originally identified from the current cell type
            contamination_candidates = list(markers_all - set(markers[celltype]))
            # identify sources of contamination
            celltype_contamination_genes = [gene for gene in contamination_candidates
                                            if pvals.loc[celltype, gene] < -np.log10(pval_cutoff) * (-1)]
            contamination_genes[celltype] = celltype_contamination_genes
            print('{}: filtered {} genes: '.format(celltype, len(celltype_contamination_genes)), flush=True)
            print(', '.join(celltype_contamination_genes), flush=True)
            print('\n', flush=True)
            # sanitize matrix for the current cell type
            celltype_idxs = self.annotations.loc[self.annotations.cell_type == celltype].index
            self.sanitized_expr.loc[celltype_idxs, celltype_contamination_genes] = 0

        if save_intermediates or (save_intermediates is None and self.save_intermediates):
            print('Saving sanitized expression matrix...')
            self.sanitized_expr.to_csv(os.path.join(self.save_dir, '{}_cell_by_gene_sanitized.csv'.format(self.sample)))
    
    def create_cytospace_input(self, sc_expression_file, sc_annotation_file, col_sc_annotation='CellType', sc_overlap_count=50, run_by_celltype=True):
        """
        Creates input files for CytoSPACE, given relevant pd.DataFrames.

        sc_expression_file  (str) : TSV filename of gene x cell scRNA-seq expression matrix, with gene names as index
        sc_annotation_file  (str) : TSV filename of scRNA-seq annotations, with cell IDs as index and a cell type column
        col_sc_annotation   (str) : Name of column in sc_annotation_file that stores the cell types
        sc_overlap_count    (int) : Will only retain cells in scRNA-seq that express >= x genes in the MERSCOPE panel.
                                    50 by default.
        run_by_celltype    (bool) : If True, will save separate input files for each cell type,
                                        so that CytoSPACE can be run separately for each cell type.
                                    If False, will save a single set of input files for a single run of CytoSPACE on the entire sample.
                                    True by default.

        Returns None; creates 5 TSV files for input to CytoSPACE.
        """
        ## scRNA-seq dataset
        sc_expressions = pd.read_csv(sc_expression_file, sep = '\t', index_col = 0)
        sc_annotations = pd.read_csv(sc_annotation_file, sep = '\t', index_col = 0)

        # retain cells whose cell types intersect with those in MERSCOPE
        common_celltypes = set(self.cell_types).intersection(sc_annotations[[col_sc_annotation]].unique())
        sc_annotations = sc_annotations.loc[sc_annotations[[col_sc_annotation]].isin(common_celltypes), :]

        # retain cells that express >= sc_overlap_count genes
        common_genes = set(self.gene_list).intersection(sc_expressions.index)
        sc_subset = sc_expressions.loc[list(common_genes)]
        sc_subset = sc_subset.loc[:, (sc_subset>0).sum(axis=0) >= sc_overlap_count]

        cells_retained = sc_annotations.index.intersection(sc_subset.columns)
        sc_expressions = sc_expressions.loc[:, cells_retained]
        sc_annotations = sc_annotations.loc[cells_retained, :]

        ## MERSCOPE dataset
        # format count matrix
        if self.sanitized_expr is not None:
            st_expressions = self.sanitized_expr.T
        else:
            print('Sanitized expression matrix does not exist; using raw expression matrix for CytoSPACE input...')
            st_expressions = self.raw_expr.T
        if st_expressions.columns.dtype == int:
            # add prefix to column names (cell IDs) if they are simple integers
            st_expressions = st_expressions.add_prefix('cell_')
        st_expressions.index.name = 'GENES'

        # create coordinate table
        st_coordinates = self.metadata
        if st_coordinates.index.dtype == int:
            # add prefix to index (cell IDs) to match st_expressions
            st_coordinates = st_coordinates.rename('cell_{}'.format)
        st_coordinates = st_coordinates.loc[st_expressions.columns, ['Y', 'X']]
        st_coordinates.columns = ['row', 'col']
        st_coordinates.index.name = 'SpotID'

        # initialize cell type fraction table
        st_annotations = self.annotations
        if st_annotations.index.dtype == int:
            # add prefix to index (cell IDs) to match st_expressions
            st_annotations.index = 'cell_' + st_annotations.index.astype(str)

        # save files
        if run_by_celltype:
            for celltype in sc_annotations[[col_sc_annotation]].unique():
                selected_cells = st_annotations.loc[st_annotations.cell_type == celltype].index
                if len(selected_cells) == 0:
                    continue
                
                # create subdirectory for this celltype
                dirname = os.path.join(self.save_dir, 'subset_{}'.format(celltype))
                Path(dirname).mkdir(exist_ok=True)

                st_expressions.loc[:, selected_cells].to_csv(os.path.join(dirname, '{}_{}'.format(self.sample, CYTOSPACE_ST_EXPRESSION_SUFFIX)), sep='\t')
                st_coordinates.loc[selected_cells, :].to_csv(os.path.join(dirname, '{}_{}'.format(self.sample, CYTOSPACE_ST_COORDINATE_SUFFIX)), sep='\t')
                celltypes = st_annotations.loc[selected_cells, :]
                celltype_list = sorted(celltypes.celltype.unique(), key=str.lower) # this will just be the current celltype
                st_annotations_subset = celltypes.celltype.value_counts().reindex(celltype_list, fill_value=0)
                st_annotations_subset = st_annotations_subset / sum(st_annotations_subset)
                st_annotations_subset = pd.DataFrame(st_annotations_subset).T
                st_annotations_subset.index = ['Fraction']
                st_annotations_subset.index.name = 'Index'
                st_annotations_subset.to_csv(os.path.join(dirname, '{}_{}'.format(self.sample, CYTOSPACE_ST_CELLTYPE_SUFFIX)), sep='\t')

                selected_sc = sc_annotations.loc[sc_annotations[[col_sc_annotation]] == celltype].index
                sc_expression.loc[:, selected_sc].to_csv(os.path.join(dirname, '{}_{}'.format(self.sample, CYTOSPACE_SC_EXPRESSION_SUFFIX)), sep='\t')
                sc_celltypes.loc[selected_sc, :].to_csv(os.path.join(dirname, '{}_{}'.format(self.sample, CYTOSPACE_SC_CELLTYPE_SUFFIX)), sep='\t')

                print('Saved CytoSPACE input files for cell type: {}'.format(celltype))
        else:
            st_annotations = st_annotations.loc[st_expressions.columns, :].cell_type.value_counts()\
                                .reindex(sorted(sc_annotations[[col_sc_annotation]].unique()), fill_value=0)
            st_annotations = st_annotations / sum(st_annotations)
            st_annotations = pd.DataFrame(st_annotations).T
            st_annotations.index = ['Fraction']
            st_annotations.index.name = 'Index'

            sc_expressions.to_csv(os.path.join(self.save_dir, '{}_{}'.format(self.sample, CYTOSPACE_SC_EXPRESSION_SUFFIX)), sep='\t')
            sc_annotations.to_csv(os.path.join(self.save_dir, '{}_{}'.format(self.sample, CYTOSPACE_SC_CELLTYPE_SUFFIX)), sep='\t')
            st_expressions.to_csv(os.path.join(self.save_dir, '{}_{}'.format(self.sample, CYTOSPACE_ST_EXPRESSION_SUFFIX)), sep='\t')
            st_coordinates.to_csv(os.path.join(self.save_dir, '{}_{}'.format(self.sample, CYTOSPACE_ST_COORDINATE_SUFFIX)), sep='\t')
            st_annotations.to_csv(os.path.join(self.save_dir, '{}_{}'.format(self.sample, CYTOSPACE_ST_CELLTYPE_SUFFIX)), sep='\t')

            print('Saved CytoSPACE input files!')
