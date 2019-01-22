""" This class sets the boundary conditions needed for a simulation run.
"""
import numpy as np
import os
from collectResults import CollectResults

class SetParam(object):
    """ ATTRIBUTES
    """
    def __init__(self):
        self.class_name = 'SetParam'
    """ METHODS
    """
    def WriteCfg(self):
        """ This method automatically generate a configuration file for the SU2 simulations
        """
        pd = list()
        pd.append( '% Configuration File for SU2 simulation:')
        pd.append( '% Part 1: PROBLEM DEFINITION'    )
        pd2 = input("Phisical problem? (EULER, NAVIER_STOKES,WAVE_EQUATION, HEAT_EQUATION, FEM_ELASTICITY, POISSON_EQUATION\n")
        pd.append( str('PHISICAL_PROBLEM='+ pd2))
        pd3 = input("Kind of Turbulent Model? (NONE, SA, SA_NEG, SST, SA_E, SA_COMP, SA_E_COMP) \n")
        pd.append( str('KIND_TURB_MODEL='+pd3))
        pd4 = input("Regime type (COMPRESSIBLE, INCOMPRESSIBLE)\n")
        pd.append( str('REGIME_TYPE='))
        pd5 = input(" Axisymmetric simulation, only compressible flows (NO, YES)\n")
        pd.append( str('AXISYMMETRIC='+pd5))
        pd6 = input("Restart solution (NO, YES) \n")
        pd.append( str('RESTART_SOL='+pd6))
        pd7 = input("Discard the data storaged in the solution and geometry files e.g. AOA, dCL/dAoA, dCD/dCL, iter, etc. (YES, NO) \n")
        pd.append( str('DISCARD_INFILES='+pd7))
        pd8 = input ("System of measurements (SI, US) \n")
        pd.append( str('SYSTEM_MEASUREMENTS='+pd8))
        pd.append( str('%----------------------------------------------------------------------------------------------%'))
        
        with open('config.cfg', 'w')  as d:
            for i in range(10):
                d.write(str('\n'+pd[i]))
                """
                d.write(pd0)
                d.write('\n'+pd1)
                d.write('\n'+pd2)
                d.write('\n'+pd3)
                d.write('\n'+pd4)
                d.write('\n'+pd5)
                d.write('\n'+pd5)
                d.write('\n'+pd6)
                d.write('\n'+pd7)
                d.write('\n'+pd8)
                d.write('\n'+pd9)
                """

    


