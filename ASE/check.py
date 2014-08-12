#!/usr/bin/python

import os
from errors_handling import *


def iffloat(var, msg):
    try:
        time_step = float(var)
    except:
        raise InputError(var,msg)

def ifstring(var, msg):
    try:
        nstep = str(var)
    except:
        raise InputError(var,msg)

def ifinteger(var,msg):
    try:
        nstep = int(var)
    except:
        raise InputError(var,msg)

def iffile(path):
    print path
    if not os.path.isfile(path):
        raise FileNotFoundError(path)

def ifdir(path):
    if not os.path.isdir(path):
        raise DirectoryNotFoundError(path)






if __name__ == '__main__':
    pass
