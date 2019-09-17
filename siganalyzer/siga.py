import sys
import argparse
import os

from .bnmf import ardnmf

def main():
    parser = argparse.ArgumentParser(description='Signature Analyzer Redux.')
    parser.add_argument('-i', '--input', required=True, help='<Required> Input matrix (maf, expression).')
    parser.add_argument('-e','--expression', help='Input is expression matrix.', action='store_true')
    parser.add_argument('-o','--output_dir', help='Output directory.', default=None)
    parser.add_argument('-c','--cosmic_signatures', help='For use with mutational signature analysis.',
                        default='../ref/cosmic_v2/cosmic_v2.txt')
    parser.add_argument('-n','--n_runs', help='Number of iterations to run ARD-NMF.', default=10)

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

    # TODO: implement wrapper

    for n_iter in args.n_runs:




if __name__ == "__main__":
    main()
