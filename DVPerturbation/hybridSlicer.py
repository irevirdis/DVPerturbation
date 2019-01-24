""" The objects have two main matrixes: a rough starting surface and a matrix which contains the connectivities  
"""
import numpy as np
import os
from meshDeform import MeshDeform
from collectResults import CollectResults

""" This class has been written to deform by applying a HicksHenne in 3D
"""
class HybridSlicer(object):
    """ ATTRIBUTES
    """
    def __init__(self, mesh=None, surface=None):
        self.class_name = 'HybridSlicer'
        self.surface = surface
        self.mesh    = mesh
        # PART 1: 
        # extract the matrix of the points inside the attribute 'surface'
        obj     = MeshDeform(mesh_name=self.mesh)
        matrix1 = obj.ExtractSurface(surface_name=self.surface)

        # check part
        with open('draft_of_hybrid.txt', 'w') as d:
            d.write(str(matrix1))

        # PART 2: 
        # section MARKER_ELEMS: c0 = connectivity, c1=node1, c2=node2, c3=node3, c4=optional|quad cells
        #print 'Initialization of HybridSlicer class: connectivity matrix (MARKER_ELEMS) is extracted from the total mesh'
        start_read  = list() 
        upper_limit = list()
        rows = [x for x in open(self.mesh).readlines()] 
        dimension=0
        limit=0
        for i in range(len(rows)):
            if 'MARKER_TAG' and self.surface in rows[i]:
                #print 'the row is:', rows[int(i+1)]
                start_read.append(int(i+2))
                limit = (int(filter(str.isdigit, rows[int(i+1)])))
                #print 'the number of cought elements is: (in python index)', limit
                k = int(i+1) 
                upper_limit.append( int(k+limit))
                #print 'the last line that will read is: (in python index)', upper_limit
                elem = rows[int(i+2)].split()
                dimension = int(len(elem)-1)
                print 'the number of the columns is:', dimension
                print 'connection is:', 
        #print 'sart:', start_read 
        #print 'upper:',upper_limit 
        print 'start_read list is:', (start_read)
        print 'upper limit list is:', (upper_limit)
        print 'limit', limit, 'dimension', dimension
        connection = np.zeros((limit, int(4*dimension)))
        for i in range(len(rows)): 
            if (i>= start_read[0]) and (i<=upper_limit[0]):
                elem = rows[i].split()
                for j in range(int (len(elem)-1)):
                    print 'the elem is:' , elem[int(j+1)]
                    index = int(i-start_read[0])
                    print 'the index is:', index
                    print 'the elem j+1 is:', elem[j+1]
                    connection[index,j] = elem[j+1]
        with open('draft_connections.txt', 'w') as d:
            d.write(str(connection))

        os.system("mkdir "+str(self.class_name+"_OUTPUTs"))
        description = 'this class is hybridslicer; this doc must be completed'
        involved_outputs=np.array(['draft_of_hybrid.txt', 'draft_connections.txt'])
        invoked_method = 'init'
        collection = CollectResults(class_name=self.class_name, involved_outputs=involved_outputs, invoked_method=invoked_method, description=description)    
        collection.Collect()


