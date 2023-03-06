import argparse

def add_plotting_arguments(parser):
    parser.add_argument("-g", "--geometry", help="ST geometry, either 'honeycomb' or 'square' accepted", type=str,
                        default="honeycomb")
    parser.add_argument("-nc", "--num-column", help="Number of columns in figure", type=int,
                        default=3)
    parser.add_argument("-mp", "--max-num-cells-plot", help="Maximum number of cells to plot in single-cell visualization", type=int, default=50000)


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

    # I/O options

    parser.add_argument("-stp", "--st-path",
                          help="Path to spatial transcriptomics data (expressions)", type=str,
                          default=None)
    parser.add_argument("-cp", "--coordinates-path",
                          help="Path to transcriptomics data (coordinates)", type=str,
                          default=None)
    parser.add_argument("-srp", "--spaceranger-path",
                          help="Path to SpaceRanger tar.gz data file", type=str,
                          default=None)

    parser.add_argument("-stctp", "--st-cell-type-path", type=str,
                          help="Path to ST cell type file (recommended for single-cell ST)",
                          default=None)
    parser.add_argument("-ctfep", "--cell-type-fraction-estimation-path", type=str,
                          help="Path to ST cell type fraction file (recommended for bulk ST)",
                          default=None)
    parser.add_argument("-ncpsp", "--n-cells-per-spot-path",
                        help="Path to number of cells per ST spot file", type=str,
                        default=None)
    parser.add_argument("-o", "--output-folder", help="Relative path to the output folder",
                        type=str,
                        default="cytospace_results")
    parser.add_argument("-op", "--output-prefix", type=str, default="",
                        help="Prefix of results stored in the 'output_folder'")

    # Solver / method options
    parser.add_argument("-mcn", "--mean-cell-numbers", type=int, default=5,
                        help="Mean number of cells per spot, default 5 (appropriate for Visium). If analyzing legacy spatial transcriptomics data, set to 20")
    parser.add_argument("-sc", "--single-cell", help="Use single-cell spatial approach if specified", action="store_true")
    parser.add_argument("-noss", "--number-of-selected-spots", help="Number of selected spots from ST data used in eahc iteration", type=int,
                        default=10000)    
    parser.add_argument("-sss", "--sampling-sub-spots", help="Sample subspots to limit the number of mapped cells if specified", action="store_true")
    parser.add_argument("-nosss", "--number-of-selected-sub-spots", help="Number of selected subspots from ST data to limit the number of mapped cells", type=int,
                        default=10000)
    parser.add_argument("-nop", "--number-of-processors", type=int, default=4,
                        help="Number of processors used for the analysis")

    parser.add_argument("-sm", "--solver-method", default="lapjv",
                        help="Which solver to use for the linear assignment problem, default 'lapjv'",
                        choices=["lapjv", "lapjv_compat", "lap_CSPR"])
    parser.add_argument("-dm", "--distance-metric", default="Pearson_correlation",
                        help="Which distance metric to use for the cost matrix, default 'Pearson_correlation'",
                        choices=["Pearson_correlation", "Spearman_correlation", "Euclidean"])
    parser.add_argument("-sam", "--sampling-method", default="duplicates",
                        help="Which underlying method to use for dealing with duplicated cells, default 'duplicates'",
                        choices=["duplicates", "place_holders"])
    parser.add_argument("-se", "--seed", help="Set seed for random generators, default 1", type=int,
                        default=1)


    # Plotting
    parser.add_argument("-p", "--plot-off", help="Turn create plots on/off", action="store_true")
    add_plotting_arguments(parser)

    arguments = parser.parse_args()

    return arguments.__dict__

