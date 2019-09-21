import sys
import argparse
import os
import pkg_resources
import pandas as pd
from .spectra import get_spectra_from_maf

from .bnmf import ardnmf

def cmd(**nmf):

def main():
    parser = argparse.ArgumentParser(description='Signature Analyzer Redux.')
    parser.add_argument('-i', '--input', required=True, help='<Required> Input matrix (maf, expression).')
    parser.add_argument('-e','--expression', help='Input is expression matrix.', action='store_true')
    parser.add_argument('-o','--output_dir', help='Output directory.', default=None)
    parser.add_argument('-c','--cosmic_signatures', help='Cosmic signatures to map to.',
                        default='cosmic2', choices=['cosmic2','cosmic3'])
    parser.add_argument('--hg', help='Hg build for mapping contexts.', default=None, choices=['hg19', 'hg38', None])
    parser.add_argument('-n','--n_runs', help='Number of iterations to run ARD-NMF.', default=10)
    parser.add_argument('-v','--verbose', help='Verbosity.', default=False, action='store_true')

    # NMF Options
    parser.add_argument('--K0', help='Initial K parameter', required=False, default=None, type=int)
    parser.add_argument('--max_iter', help='maximum iterations', required=False, default=10000, type=int)
    parser.add_argument('--del_', help='Early stop condition based on lambda change', required=False, default=1,
                        type=int)
    parser.add_argument('--tolerance', help='Early stop condition based on max lambda entry', required=False, default=1e-6,
                        type=float)
    parser.add_argument('--phi', help='dispersion parameter see paper for discussion of choosing phi '
                                      'default = 1', required=False, default=1.0, type=float)
    parser.add_argument('--a', help='Hyperparamter for lambda. We recommend trying various values of a. Smaller values'
                                    'will result in sparser results a good starting point might be'
                                    'a = log(F+N)', required=False, default=10.0, type=float)
    parser.add_argument('--b', help='Hyperparamter for lambda. Default used is as recommended in Tan and Fevotte 2012',
                        required = False,type=float, default = None)
    parser.add_argument('--objective',help='Defines the data objective. Choose between "poisson" or "gaussian". Defaults to Poisson',
                        required=False,default='poisson',type=str)
    parser.add_argument('--prior_on_W',help = 'Prior on W matrix "L1" (exponential) or "L2" (half-normal)'
                        ,required = False, default = 'L1',type=str)
    parser.add_argument('--prior_on_H',help = 'Prior on H matrix "L1" (exponential) or "L2" (half-normal)'
                        ,required = False, default = 'L1',type=str)
    parser.add_argument('--labeled', help='Input has row and column labels', required=False,default=False, action='store_true')
    parser.add_argument('--report_frequency', help='Number of iterations between progress reports', required=False,
                        default=100, type=int)

    args = parser.parse_args()

    # -------------------------------------
    # Load inputs
    # -------------------------------------
    if args.expression:
        raise Exception("Not yet implemented.")

    else:
        os.makedirs(args.output_dir, exist_ok=True)

        if args.hgfile == 'hg19':
            args.hgfile = pkg_resources.resource_filename('siganalyzer', './ref/hg19.2bit')
        else:
            args.hgfile = pkg_resources.resource_filename('siganalyzer', './ref/hg38.2bit')

        if args.cosmic_signatures == 'cosmic2':
            cosmic = pd.read_csv(pkg_resources.resource_filename('siganalyzer', './ref/cosmic_v2/cosmic_v2.txt'), sep='\t').dropna(1)
        elif args.cosmic_signatures == 'cosmic3':
            raise Exception("Not yet implemented. Provides support for only cosmic2.")
        else:
            raise Exception("Not yet implemented. Provides support for only cosmic2.")

        maf,spectra = get_spectra_from_maf(pd.read_csv(args.input, sep='\t'), hgfile=args.hgfile)
        cmap = maf[["context96.num","context96.word"]].set_index("context96.num").drop_duplicates()




        for n_iter in range(args.n_runs):








            ardnmf()





if __name__ == "__main__":
    main()
