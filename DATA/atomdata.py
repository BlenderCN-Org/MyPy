from MyPy.ERRORS.errors_handling import *

class atomdata():

    def __init__(self):
        self.weight = {
            'c': 12.01,
            'n': 14.01,
            'o': 16.00,
            'h':  1.01,
            's': 32.06,
        }
        
        self.atnum = {
            'c':  6,
            'n':  7,
            'o':  8,
            'h':  1,
            's': 16,
        }

        self.atnum_r = {}
        
        for a in self.atnum.items():
            self.atnum_r[a[1]] = a[0]

    def get_weight(self,atom):
        try: 
            return self.weight[self.atnum_r[int(atom)]]
        except:
            return self.weight[atom.lower()]

    def transform_at_numb(self,atom):
        try:
            return self.atnum_r[atom]
        except:
            raise ImplementationError(atom,'Not found in atomdata.py')
            
            
            
