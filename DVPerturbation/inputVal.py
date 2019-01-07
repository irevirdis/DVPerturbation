""" This class elaborates the inputs of the user.
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
"""
class InputVal(object):
	""" ATTRIBUTES
	"""
	def __init__(self, bump=None, positions=None):
		self.B = np.matrix(bump)
		self.pos = np.matrix(positions)
		# the following attribute will be the output given by the class "set_objective_constraints" in the future.
		# for the this first draft it is just a local string variable; constraints are absent.
		self.objective = 'Drag Coefficient'
		rows    = len(self.B)
		columns = int(self.B.size)/rows 
		if (len(self.pos) != columns):
			print '                                          WARNING !'
			print '-----------------------------------------------------------------------------------------------------------------------------'
			print 'There is an error in the definition of the input:\nThe number of control points given by B is not the same of positions file'
			print '-----------------------------------------------------------------------------------------------------------------------------'
		self.simulations = rows
		self.control_pnts = columns
	
	""" METHODS
	"""
	def MethodDiscrimination(self):
		""" Depending on the input given by the user, it will be set the 
		    numerical resolution of the gradients.
		    HYPOTHESIS: 
			1) an Identity Bump matrix stands for a step-by-step calculation of the gradients.
			2) a matrix R different from an Identity stands for a verification of the gradient direction.
		"""
		rows    = len(self.B)
		columns = int(self.B.size)/rows
		print 'the number of configurations is:', columns
		print 'the simulations will be repeated', rows, 'times'
		checking_matrix = np.ones((rows, columns))
		#print 'the checking matrix is:', checking_matrix
		method = np.count_nonzero(np.multiply(checking_matrix, self.B)) == 0
		return method
	
	def SetSimNumber(self):
		obj = CfgGenerator(self.simulations, self.B, self.pos)
		obj.WriteDraft()
		print 'from input clas: the number of simulations should be in the previous line!'
		
