#!/usr/bin/python

import sys

class percentage_with_fan():

    def __init__(self,total):
        self.total = total
        self.bar(0)
        

    def bar(self,partial):
        percentage = float(partial)/float(self.total)*100
        advancement = int(percentage)
        fan = ['-', '\\', '|', '/']
        sys.stderr.write("\r[%-100s]%3.1f%% %1s" % tuple([advancement*'=',percentage,fan[int(partial)%4]]))
        sys.stderr.flush()
