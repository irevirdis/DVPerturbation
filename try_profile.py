""" The aim of the present class is to deform the mesh nodes, starting from an su2 format.
"""
import numpy as np
import os
import math
from decimal import Decimal
from DVPerturbation import *


#filename='profile.txt'
#row = [x for x in open(filename).readlines()]
#points = np.zeros((len(row),2))
#print points.shape

md=MeshDeform('mesh_RAE2822_turb.su2')
pp=md.ExtractSurface('AIRFOIL')
print pp
print md.indices
print md.points[:,0:2]
points=md.points[:,0:2]

"""
for j in range(len(row)):
    elem=np.array(row[j].split())
    points[j,0]=np.float(elem[0])
    points[j,1]=np.float(elem[1])
"""

obj = Profile(points)

b = np.matrix([[1.00, 0.0100],
               [0.20, 0.0000],
               [0.30, 0.00],
               [0.40, 0.0000],
               [0.80, 0.0000]])
obj.ApplyBump(bump_array= b)
obj.ExportDeformed()

