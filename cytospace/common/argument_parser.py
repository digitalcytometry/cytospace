import argparse

def add_plotting_arguments(parser):
    parser.add_argument("-nr", "--num-row", help="Number of rows in pdf figure", type=int,
                        default=4)
    parser.add_argument("-nc", "--num-column", help="Number of columns in pdf figure", type=int,
                        default=4)
    parser.add_argument("-r", "--rotation-flag", help="Rotate plot", action="store_false")
    parser.add_argument("-nv", "--plot-nonvisium", help="Plot with custom slide dimensions", action="store_true")
    parser.add_argument("-rd", "--rotation-degrees", help="Rotation on plot", type=int,
                        default=270)
    parser.add_argument("-ss", "--spot-size", help="Set size of ST spots", type=int,
                        default=175)
    parser.add_argument("-pm", "--plot-marker", help="Shape of ST spots", type=str,
                        default='h')


def argument_parser():
    parser = argparse.ArgumentParser(description="CytoSPACE is a computational strategy for "
                                                 "assigning single-cell transcriptomes to in situ "
                                                 "spatial transcriptomics (ST) data. Our "
                                                 "method solves single cell/spot assignment by "
                                                 "minimizing a correlation-based cost function "
                                                 "through a linear programming-based optimization "
                                                 "routine.")

    # Required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument("-sp", "--scRNA-path", help="Path to scRNA-Seq data", type=str,
                          default=None, required=True)
    required.add_argument("-ctp", "--cell-type-path", help="Path to cell type labels", type=str,
                          default=None, required=True)
    required.add_argument("-stp", "--st-path",
                          help="Path to spatial transcriptomics data (expressions)", type=str,
                          default=None, required=True)
    required.add_argument("-cp", "--coordinates-path",
                          help="Path to transcriptomics data (coordinates)", type=str,
                          default=None, required=True)
    required.add_argument("-ctfep", "--cell-type-fraction-estimation-path", type=str,
                          help="Path to cell type fraction file",
                          default=None, required=True)

    # I/O options
    parser.add_argument("-ncpsp", "--n-cells-per-spot-path",
                        help="Path to number of cells per ST spot file", type=str,
                        default=None)
    parser.add_argument("-o", "--output-folder", help="Relative path to the output folder",
                        type=str,
                        default="cytospace_results")
    parser.add_argument("-op", "--output-prefix", type=str, default="",
                        help="Prefix of results stored in the 'output_folder'")
    parser.add_argument("-d", "--delimiter", type=str, default="\t",
                        help="Set delimiter of the input files, default '\\t'")

    # Solver / method options
    parser.add_argument("-sm", "--solver-method", default="lapjv",
                        help="Which solver to use for the linear assignment problem, default 'lapjv'",
                        choices=["lapjv", "lapjv_compat", "lap_CSPR"])
    parser.add_argument("-sam", "--sampling-method", default="duplicates",
                        help="Which unerlying method to use for dealing with duplicated cells, default 'duplicates'",
                        choices=["duplicates", "place_holders"])
    parser.add_argument("-mcn", "--mean-cell-numbers", type=int, default=5,
                        help="Mean number of cells per spot, default 5 (appropriate for Visium). If analyzing legacy spatial transcriptomics data, set to 20")
    parser.add_argument("-se", "--seed", help="Set seed for random generators, default 1", type=int,
                        default=1)
    parser.add_argument("-dm", "--distance-metric", default="Pearson_correlation",
                        help="Which distance metric to use for the cost matrix, default 'Pearson_correlation'",
                        choices=["Pearson_correlation", "Spearman_correlation", "Cosine", "Euclidean"])
    parser.add_argument("-nosc", "--number-of-selected-cells", help="Number of selected cells from scRNA-seq data used in eahc iteration", type=int,
                        default=10000)
    parser.add_argument("-noss", "--number-of-selected-spots", help="Number of selected spots from ST data used in eahc iteration", type=int,
                        default=10000)
    parser.add_argument("-nop", "--number-of-processors", type=int, default=4,
                        help="Number of processors used for the analysis")
    parser.add_argument("-sc", "--single-cell", help="Use single-cell spatial approach or not", action="store_true")

    # Plotting
    parser.add_argument("-p", "--plot-off", help="Turn create plots on/off", action="store_true")
    add_plotting_arguments(parser)

    arguments = parser.parse_args()

    return arguments.__dict__

