"""
process_merscope.py

Master processing script for MERSCOPE datasets.
"""

from merscope_expression import MerscopeExpression

import argparse

def process_merscope(sample, data_dir, annotation_file,
                        sc_expression_file, sc_annotation_file,
                        save_dir=None, sanitize_matrix=True, save_intermediates=True):
    obj = MerscopeExpression(sample, data_dir, save_dir, save_intermediates)
    obj.read_raw_data()
    obj.read_metadata()
    obj.set_annotations(annotation_file, 'celltype', is_abs_path=True) # assumes annotation_file is a CSV with column "celltype"
    if sanitize_matrix:
        obj.get_zplanes([0, 3])
        obj.sanitize_matrix(z_peri_idx=0, z_mid_idx=3)
    obj.create_cytospace_input(sc_expression_file, sc_annotation_file)

def main():
    parser = argparse.ArgumentParser(description='Process MERSCOPE data')
    parser.add_argument('sample', help='sample name')
    parser.add_argument('data_dir', help='path to the root directory for results, containing cell_by_gene.csv and other files')
    parser.add_argument('annotation_file', help='CSV file containing cell type annotations for the MERSCOPE data')
    parser.add_argument('sc_expression_file', help='TSV file of gene expression matrix for reference single-cell dataset')
    parser.add_argument('sc_annotation_file', help='TSV file of cell type annotations for reference single-cell dataset')
    parser.add_argument('--save-dir', help='path to save outputs, by default, saves under the current directory')
    parser.add_argument('--no-sanitization', action='store_false', help='does not sanitize expression matrix if flag is specified')
    parser.add_argument('--save-intermediates', action='store_true', help='stores all intermediate outputs if flag is specified')
    args = parser.parse_args()

    process_merscope(args.sample, args.data_dir, args.annotation_file,
                        args.sc_expression_file, args.sc_annotation_file,
                        args.save_dir, args.no_sanitization, args.save_intermediates)

if __name__ == '__main__':
    main()
