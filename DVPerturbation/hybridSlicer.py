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
        with open('draft_of_hybrid.txt', 'r') as r:
            with open('coordinates_matrix.txt', 'w') as w:
                data = r.read()
                data = data.replace("[","")
                data = data.replace("]","")
                data = data.replace("'", "")
                w.write(data)
        #os.system("rm draft_of_hybrid.txt")

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

        dim_m = len(matrix1)
        dim_c = int(upper_limit[0] - start_read[0])
        skip  = int(dim_c - dim_m)
        connection = np.zeros((limit, int(4*dimension)))

        for i in range(len(rows)):
            #m_rows=len(matrix1)
            if (i>= start_read[0]) and (i< (int(upper_limit[0]))):
                elem = rows[i].split()
                for j in range(int (len(elem)-1)):
                    #print 'the elem is:' , elem[int(2*(j+1))]
                    index = int(i-start_read[0])
                    #print 'the index is:', index
                    #print 'the elem j+1 is:', elem[j+1]
                    connection[index,j] = elem[j+1]
                    # part dedicated to the x, y, z  coordiantes:
                    linked_row = connection[index,j]
                    c_index = int(3*j + len(elem)-1)
                    c_index_x = c_index
                    c_index_y = int(c_index + 1)
                    c_index_z = int(c_index + 2)
                    
                    for k in range(len(matrix1)):
                        if connection[index,j] == float(matrix1[k,0]):
                            connection[index, c_index_x] = matrix1[k,1]
                            connection[index, c_index_y] = matrix1[k,2]
                            connection[index, c_index_z] = matrix1[k,3]
       
        with open('draft_connections.txt', 'w') as d:
            np.savetxt('draft_connections.txt', connection, fmt='%.6e')
            #d.write(str(connection))
        with open('draft_connections.txt', 'r') as r:
            with open('connections_matrix.txt', 'w') as w:
                data = r.read()
                data = data.replace("[","")
                data = data.replace("]","")
                data = data.replace("'", "")
                w.write(data)
        #os.system("rm draft_connections.txt")

        os.system("mkdir "+str(self.class_name+"_OUTPUTs"))
        description = 'this class is hybridslicer; this doc must be completed'
        involved_outputs=np.array(['coordinates_matrix.txt', 'connections_matrix.txt', 'draft_of_hybrid.txt', 'draft_connections.txt'])
        invoked_method = 'init'
        collection = CollectResults(class_name=self.class_name, involved_outputs=involved_outputs, invoked_method=invoked_method, description=description)    
        collection.Collect()


