#!/usr/bin/python

class FileNotFoundError(Exception):
    """Exception raised for missing files.

    Attributes:
        message --- name of the missing file
    """
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)

class DirectoryNotFoundError(Exception):
    """Exception raised for missing directory.

    Attributes:
        message --- name of the missing file
    """
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)


class InputError(Exception):
    """Exception raised when there is an input error."
    
    Attributes:
        expression --- input expression in which the error occurred
        message --- explanation of the error
    """
    
    def __init__(self,expression,message):
        self.expression = expression
        self.message = message
    def __str__(self):
        return repr('%s -> %s' % tuple([self.expression,self.message]))


class ImplementationError(InputError):
    """Ecveption raised when something is defined only partially inside the code.

    Attributes:
       expression --- requested espression
       message --- what is the problem
    """
