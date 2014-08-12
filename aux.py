#!/bin/bash


from errors_handling import *


try:
    time_step = float(time_step)
except:
    raise InputError(time_step,'The time step has to be a float!')







if __name__ == '__main__':
    pass
