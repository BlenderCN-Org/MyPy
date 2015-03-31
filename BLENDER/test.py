#!/usr/bin/python

import pybel
import json
# import sys

# def main():
#     # Use openbabel to read the molecule format
#     molecule = pybel.readfile('xyz',sys.argv[1])


def molecule_to_json(molecule):
    """Converts an OpenBabel molecule to json for use in Blender."""
    import openbabel

    # Save atom element type and 3D location.
    atoms = [{"element": atom.type,
              "location": atom.coords}
             for atom in molecule.atoms]

    # Save number of bonds and indices of endpoint atoms
    bonds = [{"atoms": [b.GetBeginAtom().GetIndex(), b.GetEndAtom().GetIndex()],
              "order": b.GetBondOrder()}
             for b in openbabel.OBMolBondIter(molecule.OBMol)]

    return json.dumps({"atoms": atoms, "bonds": bonds}, indent=4)



# if __name__ == '__main__':
#     main()




molecule = pybel.readstring("smi", "O=C1C2=C(N=CN2C)N(C(=O)N1C)C")
molecule.make3D()

with open("caffeine.json", "w") as out_file:
    out_file.write(molecule_to_json(molecule))
