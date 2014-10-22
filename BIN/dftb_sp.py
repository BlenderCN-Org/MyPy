#!/usr/bin/python
"""Using the ASEInterface to run a geometry optimization.

   This script works directly from command line, so use the
   --help option to see the help."""

from MyPy.ERRORS.errors_handling import *


def main():
    from MyPy.ASE.ASEInterface import ASEInterface
    args = __parser()

    ase = ASEInterface(job=args.job)
    ase.setMachine(__setHostname())
    ase.setGeneralParams(args.method, args.xyzfile, charge = args.charge,
                         spinm = args.spinm)
    ase.setGeomOptParams('out')

    ase.getParams()
    print ase.runJob()



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
                        dest='method',
                        type=str,
                        default='dftb_std',
                        metavar = 'METHOD',
                        choices=['dftb_std', 'dftb_std-D3', 'dftb_std-dDMC',
                                 'dftb_std-D3H4', 'dftb_mio11'],
                        help='select the method to compute the forces (default: dftb_std)')

    parser.add_argument('-c', '--charge',
                        action='store', # optional because action defaults to 'store'
                        dest='charge',
                        type=int,
                        metavar='CHARGE',
                        default=0,
                        help='define system charge (default: 0)')

    parser.add_argument('-s', '--spinm',
                        action = 'store',
                        dest = 'spinm',
                        type = int,
                        metavar = 'SPINMULT',
                        default = 1,
                        help='define the spin multiplicity (default: 1)')

    parser.add_argument('xyzfile',
                        action = 'store',
                        type = str,
                        metavar = '<XYZFILE>',
                        help = 'file containing the starting structure')

    parser.add_argument('-j', '--job',
                        action = 'store',
                        dest = 'job',
                        type = str,
                        metavar = 'JOB',
                        default = 'singlePoint',
                        choices = ['singlePoint', 'geomOpt']
                        )

    return parser.parse_args()
    # (options, args) = parser.parse_args(sys.argv[:])

    # if len(args) != 0:
    #     parser.error("wrong number of arguments")

    # return options






if __name__ == '__main__':
    main()
