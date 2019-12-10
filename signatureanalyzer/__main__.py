"""
For argument parsing and CLI interface
"""
import sys
import argparse
import os
import pkg_resources
from argparse import RawTextHelpFormatter

from .signatureanalyzer import run_maf
from .signatureanalyzer import run_spectra
from .signatureanalyzer import run_matrix

def main():
    parser = argparse.ArgumentParser(description='Signature Analyzer GPU.', formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        'input',
        help="Input matrix for decomposition. Signature Analyzer uses the format of (samples x features)\n"
             "Assumes input is a .maf by default and will compute the 96-base context spectra if not provided\n"
             "  * Use {-type} to specific different input types"
    )
    parser.add_argument(
        '-t','--type',
        default='maf',
        help="Input type. Specify whether input is a .maf, a 96 base context spectra, or an RNA expression matrix (default: 'maf')\n"
             "  * NOTE: for expression is is reccomended to use log-transformed & gaussian {--objective} function",
        choices=['maf','spectra','matrix']
    )
    parser.add_argument(
        '-n','--nruns',
        help="Number of iterations to run ARD-NMF. Significant speed up if GPU is available (default: 10)",
        default=10,
        type=int
    )
    parser.add_argument(
        '-o','--outdir',
        help="Directory to save outputs (default: '.')",
        default="."
    )
    parser.add_argument(
        '--cosmic',
        help="Cosmic signatures to map to and provide results for. Support for Cosmic 2 & 3 (default: 'cosmic2')\n"
             "  * Reference: https://cancer.sanger.ac.uk/cosmic/signatures",
        default='cosmic2',
        choices=[
            'cosmic2',
            'cosmic3',
            'cosmic3_exome',
            'cosmic3_DBS',
            'cosmic3_ID',
            'cosmic3_TSB',
            ]
    )
    parser.add_argument(
        '--hg_build',
        help="Path to 2bit, human genome build for mapping mutational contexts. Required if mutational context is not provided (default: None)",
        default=None
    )
    parser.add_argument(
        '--cuda_int',
        help="GPU to use. Defaults to (cuda: 0). If (None), or if no GPU is available, will default to CPU (default: 0)",
        default=0
    )
    parser.add_argument(
        '--verbose',
        help="Verbosity",
        default=False,
        action='store_true'
    )

    # -----------------------------------------
    # NMF Options for Argparse
    # -----------------------------------------
    parser.add_argument(
        '--K0',
        help="Initial K0 parameter. If not provided, ARD-NMF starts with K0 = no. features",
        default=None,
        type=int
    )
    parser.add_argument(
        '--max_iter',
        help="Maximum number of iterations for CAVI algorithm if not reached {--tolerance} (default: 10000)",
        default=10000,
        type=int
    )
    parser.add_argument(
        '--del_',
        help="Early stop condition based on lambda change (default: 1)",
        default=1,
        type=int
    )
    parser.add_argument(
        '--tolerance',
        help="Early stop condition based on max lambda entry (default: 1e-10)",
        default=1e-10,
        type=float
    )
    parser.add_argument(
        '--phi',
        help="Dispersion parameter for CAVI. See paper for details on selection (default: 1).\n"
             "   * NOTE: If using gaussian {--objective}, scaled by the variance of input matrix",
        default=1.0,
        type=float
    )
    parser.add_argument(
        '--a',
        help="Hyperparamter for lambda. We recommend trying various values of a. Smaller values will result in\n"
             "sparser results. Reccommended starting hyperparameter: a = log(F+N). (default: 10.0)",
         default=10.0,
         type=float
     )
    parser.add_argument(
        '--b',
        help="Hyperparamter for lambda. Default is computed automatically as speicified by Tan and Fevotte 2013",
        default = None,
        type=float
    )
    parser.add_argument(
        '--objective',
        help="Objective function for ARD-NMF. (default: 'poisson')\n"
             "  * mutational signatures --> poisson (DEFAULT)\n"
             "  * log-norm expression   --> gaussian",
        default="poisson",
        choices=["poisson","gaussian"],
        type=str
    )
    parser.add_argument(
        '--prior_on_W',
        help="Prior on W matrix L1 (exponential) or L2 (half-normal) (default: 'L1')",
        default="L1",
        choices=["L1","L2"],
        type=str
    )
    parser.add_argument(
        '--prior_on_H',
        help="Prior on H matrix L1 (exponential) or L2 (half-normal) (default: 'L1')",
        default="L1",
        choices=["L1","L2"],
        type=str
    )
    parser.add_argument(
        '--report_freq',
        help="Number of iterations between progress reports (default: 250)",
        default=250,
        type=int
    )
    parser.add_argument(
        '--active_thresh',
        help="Active threshold for consdiering a threshold relevant (default: 0.01)",
        default=0.01,
        type=float
    )

    # -----------------------------------------
    # NMF Post-processing Arguments
    # -----------------------------------------
    parser.add_argument(
        '--cut_norm',
        help="Min normalized value for mean signature. Used in marker selection during post-processing (matrix). (default: 0.5)",
        default=0.5,
        type=float
    )
    parser.add_argument(
        '--cut_diff',
        help="Difference between mean selected signature and mean unselected signatures for marker selection (matrix). (default: 1.0)",
        default=1.0,
        type=float
    )

    args = parser.parse_args()
    if args.cuda_int == 'None':
        args.cuda_int = None

    print("---------------------------------------------------------")
    print("---------- S I G N A T U R E  A N A L Y Z E R  ----------")
    print("---------------------------------------------------------")

    if args.type == 'maf':
        run_maf(
            args.input,
            **vars(args)
        )
    elif args.type == 'spectra':
        run_spectra(
            args.input,
            **vars(args)
        )
    elif args.type == 'matrix':
        run_matrix(
            args.input,
            **vars(args)
        )

if __name__ == "__main__":
    main()
