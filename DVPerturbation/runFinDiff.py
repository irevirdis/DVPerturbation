""" This class has got a method to automatically run a finite difference calculation.
"""
import numpy as np
import os
from cfgGenerator import CfgGenerator

""" This class takes as inputs:
    : param B         =   matrix [NxM]N is the number of the different configurations that will be studied;
                          M is the number of the control points inside the fluid domain
    : param pos       =   matrix [Mx2] where the first column contains the integer [0] for lower side OR [1] for
                          upper side of the surface; M is the numerosity of the control points; the second column 
                          contains the non-dimensional position of the control points.
    : param who       =   it represents the choice between SU2 official scripts and IVMC scripts; the default 
                          option is set to IVMC.
"""
class RunFinDiff(object):
    """ ATTRIBUTES
    """
    def __init__(self, bump =None, positions=None, who=None):
        self.B   = np.matrix(bump)
        self.pos = np.matrix(positions)
        if (who is None) or (who is 'IVMC'):
            self.who = 'IVMC'
        else:
            self.who = who

    """ METHODS
    """ 
    def Calc(self):
        """ automatic run of DEF, CFD for the finite differences case.
        """
        if self.who is 'IVMC':      
                obj = CfgGenerator(len(self.B), self.B, self.pos)
                obj.WriteDraft()
                obj.RunSU2()
                c = obj.ReadRes()
        else:
                c = 0.
                print 'from class RunFinDiff, method Calc: to be completed'

        return c

