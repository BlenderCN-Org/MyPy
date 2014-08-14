#!/usr/bin/python
"""Using the ASEInterface to run a geometry optimization.

   This script works directly from command line, so use the 
   --help option to see the help."""

from MyPy.ERRORS.errors_handling import *


def main():
    import sys
    from MyPy.ASE.ASEInterface import ASEInterface
    args = __parser()
    print args

    ase = ASEInterface(job='geomOpt')
    ase.setMachine(__setHostname())
    ase.setGeneralParams(args.method, args.xyzfile, charge = args.charge, 
                         spinm = args.spinm)
    ase.setGeomOptParams(args.out_prefix, fmax = args.fmax)

    outf = open(args.out_prefix+'.out','w',1)
    outf.write(ase.runJob())



def __setHostname():
    import socket
    hostname = socket.gethostname()
    
    if hostname == 'icmbpriv20':
        return 'workstation'
    else:
        raise ImplementationError(hostname,'Hostname not defined in __setHostname')


def __parser():
    import argparse
    parser = argparse.ArgumentParser(version='%prog 1.0', 
                                     description='Perform geometry optimizations.')
    parser.add_argument('-m', '--method',
                        action = 'store',
                        nargs = 1,
                        dest='method',
                        type=str,
                        default='dftb_std',
                        metavar = 'METHOD',
                        choices=['dftb_std','dftb_std-D3','dftb_std-dDMC','dftb_std-D3H4'],
                        help='select the method to compute the forces (default: dftb_std)')

    parser.add_argument('-c', '--charge',
                        action='store', # optional because action defaults to 'store'
                        nargs = 1,
                        dest='charge',
                        type=int,
                        metavar='CHARGE',
                        default=0,
                        help='define system charge (default: 0)')
                        
    parser.add_argument('-s', '--spinm',
                        action = 'store',
                        nargs = 1,
                        dest = 'spinm',
                        type = int,
                        metavar = 'SPINMULT',
                        default = 1,
                        help='define the spin multiplicity (default: 1)')
    
    parser.add_argument('-o','--out-prefix',
                        action = 'store',
                        nargs = 1,
                        dest = 'out_prefix',
                        type = str,
                        metavar = '<STRING>',
                        default = 'geom_opt_',
                        help = 'prefix of the output files (default: geom_opt_)')

    parser.add_argument('-t','--fmax',
                        action = 'store',
                        nargs = 1,
                        dest = 'fmax',
                        type = float,
                        metavar = 'FMAX',
                        default = 0.001,
                        help = 'the highest value acceptable for the highest force component (default: 0.001)')

    parser.add_argument('xyzfile',
                        action = 'store',
#                        nargs = 1,
#                        dest = 'xyzfile',
                        type = str,
                        metavar = '<XYZFILE>',
                        help = 'file containing the starting structure')



    return parser.parse_args()
    # (options, args) = parser.parse_args(sys.argv[:])

    # if len(args) != 0:
    #     parser.error("wrong number of arguments")

    # return options






if __name__ == '__main__':
    main()
