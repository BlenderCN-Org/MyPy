#!/usr/bin/python2.7

from MyPy.DATA.atomdata import atomdata
from MyPy.MD.xyz_trajsplit import mdxyz

def main():
    args = __parser()
    xyz = mdxyz(args.in_xyzfile)
    data = atomdata()

    for i,t in enumerate(xyz.at_type):
        xyz.at_type[i] = data.transform_at_numb(int(t))

    xyz.write_xyz('1-'+str(len(xyz.at_type)),args.out_xyzfile)


def __parser():
    import argparse
    parser = argparse.ArgumentParser(version='%prog 1.0', 
                                     description='Perform geometry optimizations.')
    parser.add_argument('in_xyzfile',
                        action = 'store',
                        type = str,
                        metavar = '<IN_XYZFILE>',
                        help = 'file containing the starting structure')

    parser.add_argument('out_xyzfile',
                        action = 'store',
                        type = str,
                        metavar = '<OUT_XYZFILE>',
                        help = 'name of the file containing the output structure')



    return parser.parse_args()




if __name__ == '__main__':
    main()
