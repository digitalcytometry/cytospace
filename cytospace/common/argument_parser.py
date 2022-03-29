import argparse


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
    parser.add_argument("-o", "--output-folder", help="Relative path to the output folder",
                        type=str,
                        default="cytospace_results")
    parser.add_argument("-op", "--output-prefix", type=str, default="",
                        help="Prefix of results stored in the 'output_folder'")
    parser.add_argument("-d", "--delimiter", type=str, default=",",
                        help="Set delimiter of the input files, for instance '\\t'")

    # Solver / method options
    parser.add_argument("-sm", "--solver-method", default="lapjv",
                        help="Which solver to use for the linear assignment problem.",
                        choices=["lapjv", "lapjv_compat", "lap_CSPR"])
    parser.add_argument("-mcn", "--mean-cell-numbers", type=int, default=5,
                        help="Mean number of cells per spot. Default is set to 5 (Visium). If analyzing legacy spatial transcriptomics data, set to 20")
    parser.add_argument("-se", "--seed", help="Set seed for random generators", type=int,
                        default=1)

    # Plotting
    parser.add_argument("-p", "--plot-off", help="Turn create plots on/off", action="store_true")
    parser.add_argument("-nr", "--num-row", help="Number of rows in pdf figure", type=int,
                        default=4)
    parser.add_argument("-nc", "--num-column", help="Number of coulmns in pdf figure", type=int,
                        default=4)
    parser.add_argument("-r", "--rotation-flag", help="Rotate plot", action="store_false")
    parser.add_argument("-rd", "--rotation-degrees", help="Rotation on plot", type=int,
                        default=270)
    parser.add_argument("-ss", "--spot-size", help="Set size of ST spots", type=int,
                        default=175)
    parser.add_argument("-pm", "--plot-marker", help="Shape of ST spots", type=str,
                        default='h')

    arguments = parser.parse_args()

    return arguments.__dict__
