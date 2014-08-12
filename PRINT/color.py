#!/usr/bin/python
class colorprinter():

    def __init__(self): 
        # Color Definition Dictionary:
        self.color_def={ \
            'ENDC'    : '\033[0m',\
            'HEADER'  : '\033[95m',\
            'OKBLUE'  : '\033[94m',\
            'OKGREEN' : '\033[92m',\
            'WARNING' : '\033[93m',\
            'FAIL'    : '\033[91m'\
        }

    def __print_with_color(self,color,msg):
        print self.color_def[color]+msg+self.color_def['ENDC']

    def __stderr_with_color(self,color,msg):
        import sys
        sys.stderr.write(self.color_def[color]+msg+self.color_def['ENDC']+'\n')

    def cprint(self,color,msg):
        if not color in self.color_def.keys():
            self.__stderr_with_color('FAIL',"Color "+color+" unknown for the following message:\n%s" % str(msg))
        else:
            self.__print_with_color(color,msg)

    def cstderr(self,color,msg):
        if not color in self.color_def.keys():
            self.__stderr_with_color('FAIL',"Color "+color+" unknown for the following message:\n%s" % str(msg))
        else:
           self.__stderr_with_color(color,msg)

    def print_colors(self):
        for c in self.color_def.keys():
            if c != 'ENDC':
                self.cprint(c,c)
