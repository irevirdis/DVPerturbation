""" This class has got a method to automatically run an Adjoint calculation.
"""
import numpy as np
import os
from cfgGenerator import CfgGenerator

""" This class takes as inputs:

    : param pos       =   matrix [Mx2] where the first column contains the integer [0] for lower side OR [1] for
                          upper side of the surface; M is the numerosity of the control points; the second column 
                          contains the non-dimensional position of the control points.
    : param who       =   it represents the choice between SU2 official scripts and IVMC scripts; the default 
                          option is set to IVMC.
"""
class RunAdjoint(object):
    """ ATTRIBUTES
    """
    def __init__(self, bump=None, positions=None, who=None):
        self.pos = np.matrix(positions)
        if (who is None) or (who is 'IVMC'):
            self.who = 'IVMC'
        else:
            self.who = who
        self.B = bump

    """ METHODS
    """ 
    def Calc(self):
        """ automatic run of DEF, CFD for the finite differences case.
        """
        c = 'global'

        if self.who is 'IVMC':      
                obj = CfgGenerator(bump= self.B, pos=self.pos)
                obj.WriteDraft()
                obj.RunADJ()
        else:
                c = 'hello from RunADJ class'
                print 'from class RunAdjoint, method Calc: to be completed'

        return c

