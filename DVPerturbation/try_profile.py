""" The aim of the present class is to deform the mesh nodes, starting from an su2 format.
"""
import numpy as np
import os
import math
from profile import Profile
from collectResults import CollectResults
from decimal import Decimal


filename='profile.txt'
row = [x for x in open(filename).readlines()]
points = np.zeros((len(row),2))
print points.shape
for j in range(len(row)):
    elem=np.array(row[j].split())
    points[j,0]=np.float(elem[0])
    points[j,1]=np.float(elem[1])

obj = Profile(points)

b = np.matrix([[0.10, 0.0000],
               [0.20, 0.0000],
               [0.30, 0.0001],
               [0.40, 0.0000]])
obj.ApplyBump(bump_array= b)


