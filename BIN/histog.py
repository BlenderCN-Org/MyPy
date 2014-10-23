#!/usr/bin/python

import sys
import numpy as np

def main():
    args = __parser()
    
    if args.infile == None:
        filec = sys.stdin.readlines()
    else:
        filec = args.infile.readlines()


    data = []
    for l in filec:
        if l[0] != '#':
            try:
                data.append(float(l))
            except:
                pass

    data = np.array(data)

    if args.range == None:
        minmax = [data.min(), data.max()]
    else:
        minmax = args.range

    if all( x <= minmax[0] or x >= minmax[1] for x in data):
        sys.stderr.write('No data in the specified range.')
        exit(1)
        

    hist, bin_edges = np.histogram(data, args.nbin, range=minmax, density=args.normed)

    if args.probab:
        hist = hist / float(len(data))

    for b,h in zip(bin_edges,hist):
        print b,h



def __parser():
    import argparse
    parser = argparse.ArgumentParser(version='%prog 1.0', 
                                     description='Compute histogram of inputted data.')
    parser.add_argument('-b', '--bin',
                        action = 'store',
                        dest='nbin',
                        type=int,
                        metavar = 'BIN',
                        required = True,
                        help='number of bin.')

    parser.add_argument('-r', '--range',
                        action='store', # optional because action defaults to 'store'
                        dest='range',
                        type=float,
                        metavar='MIN MAX',
                        nargs = 2,
                        help='lower and upper range of the bins. Values outside are ignored.')
                        
    parser.add_argument('-n', '--normed',
                        action = 'store_true',
                        dest = 'normed',
                        default = False,
                        help='the integral over the range is 1.')

    parser.add_argument('-p', '--probability',
                        action = 'store_true',
                        dest = 'probab',
                        default = False,
                        help='the maximum value is 1.')
    
    parser.add_argument('-f', '--infile',
                        action = 'store',
                        type = file,
                        metavar = '<INPUT>',
                        required = False,
                        help = 'file containing a list of values to histogram. You can supply the values even from stdin')



    return parser.parse_args()




if __name__ == '__main__':
    main()
