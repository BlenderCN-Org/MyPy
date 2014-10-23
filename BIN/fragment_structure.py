#!/usr/bin/python2.7


def main():
    import MyPy.STRUCTURE.fragments as frg
    import numpy as np
    
    args = __parser()

    atomnos, atomcoords, cstring = frg.readxyz(args.xyzfile)

    voisins = frg.neighbours(atomnos, atomcoords)
    monomers = frg.molsplit(voisins)
    
    i = 0
    for m in monomers.keys():
        monomer = []
        for atom in monomers[m]:
            monomer.append([atomnos[atom]]+atomcoords[atom])
        fname = '%s_mono%i.xyz'%(args.out_prefix, m+1)
        #        writexyz(fname, monomer, cstring)
        print 'Mon %i molecular mass: %12.6f' % tuple([i,frg.molecular_mass(monomers[i],atomnos)])
        i += 1

    # print "Mon 1 molecular mass: ",
    # print "Mon 2 molecular mass: ",frg.molecular_mass(monomers[1],atomnos)

    # print "Mon 1: ", sorted(monomers[0])
    # print "Mon 2: ", sorted(monomers[1])

    for a in monomers.keys():
        outfn = args.out_prefix+str(a)+'.xyz'
        xyz_w = '%5d\n\n' % len(monomers[a])
        for b in sorted(monomers[a]):
            xyz_w += '%3d  %12.6f  %12.6f  %12.6f\n' % tuple([atomnos[b]]+atomcoords[b])
            
        open(outfn,'w').write(xyz_w)



    
#def __writexyz()
    


def __parser():
    import argparse
    parser = argparse.ArgumentParser(version='%prog 1.0', 
                                     description='Separe two molecules from the same xyz file.')
    
    parser.add_argument('-o','--out-prefix',
                        action = 'store',
                        dest = 'out_prefix',
                        type = str,
                        metavar = '<STRING>',
                        default = 'frag_',
                        help = 'prefix of the output files (default: geom_opt)')

    parser.add_argument('xyzfile',
                        action = 'store',
                        type = str,
                        metavar = '<XYZFILE>',
                        help = 'file containing the starting structure')



    return parser.parse_args()



if __name__ == '__main__':
    main()
