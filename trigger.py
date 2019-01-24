from DVPerturbation import *

#obj = CfgGenerator()
#obj.WorkingFolder()
#obj.WriteDraft()

# check of InputVal class
"""
B = np.eye(6)*1E-03
z = np.zeros([1,6])
print '-------------------------------------------------------------------------'
print 'from trigger file:'
print z
B = np.vstack((z,B))
print 'the declared Bump matrix is:'
print B

print 'the position array is :'
pos = [[0.0, 0.25], [0.0, 0.5], [0.0, 0.75], [1.0, 0.25], [1.0, 0.5], [1.0, 0.75]]
print pos
print '-------------------------------------------------------------------------'
"""
# instance to inputval

"""
obj2 = InputVal(bump = B, positions=pos)
print obj2.MethodDiscrimination()

obj2.SetSimNumber()
"""

# instance to cfgGenerator

"""
obj = CfgGenerator( len(B) , B, pos)
obj.RunSU2()
c = obj.ReadRes()
print c
"""

# instance to RunFinDiff --------------------------> everything is ok! 

#obj = RunFinDiff(bump=B, positions=pos, who='IVMC')
#obj.Calc()

# instance to RunADJ
"""
obj = RunAdjoint(positions=pos, who='IVMC', bump=B)
obj.Calc()
"""


# test of the first method inside meshDeform
##obj = MeshDeform('mesh_stator.su2') # instance ok
#obj = MeshDeform('mesh_RAE2822_turb.su2')
#matrix = obj.ExtractSurface('AIRFOIL') # ok

#obj = MeshDeform('mesh_ONERAM6_turb_tets.su2')
#matrix = obj.ExtractSurface('wing')

#obj.TransformMesh('AIRFOIL', translate = [0.2, 0.02])
#matrix = obj.SortedBlade('AIRFOIL')
#print matrix
#b = [[3,9],[0.001, 0.001]]
#obj.ApplyBump(surface='AIRFOIL', bump_array= b)

#b = np.matrix([[0.10, 0.0000],
#               [0.20, 0.0000],
#               [0.30, 0.0001],
#               [0.40, 0.0000]])
#obj.ApplyBump(surface='AIRFOIL', bump_array= b)

# test of the DEF application:
#obj.ApplyDEF(config_file='baseline.cfg')
#obj.VerifyIntegrity(mesh_old='mesh_RAE2822_turb.su2', mesh_new='mesh_out.su2', bump_matrix=b)


# part to test the PostPro class:
#obj = PostPro()
#obj.Gradient(file_name='surface_sens.dat')

#obj = SetParam()
#obj.WriteCfg()


# test of hybrid mesh
obj = HybridSlicer(mesh='mesh_ONERAM6_turb_tets.su2', surface='wing')






