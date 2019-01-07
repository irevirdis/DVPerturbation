""" The aim of the present class is to deform the mesh nodes, starting from an su2 format.
"""
import numpy as np
import os
import math

""" This class takes as inputs:
    : param mesh_name = it is a string that will points the mesh file.

    : RETURNS: a new mesh file in which the nodes coordinates have been changed.
"""
class MeshDeform(object):
    # attributes of the class
    
    def __init__(self, mesh_name):
        self.mesh_name = mesh_name
    
    # methods
    def ExtractSurface(self, surface_name):
        """ this method has been written to extract the coordinates of a surface,
            given the name as input.

            MAKE SURE THAT THE SURFACE NAME WILL BE WRITTEN IN LOWER OR UPPER CASE
            AS REPORTED INTO THE ORIGINAL MESH FILE.
        """
        print '--------------------------------------------------------------'
        print '         Start calculation from ExtractSurface Method'
        filename = self.mesh_name
        surface = surface_name

        row = [x for x in open(filename).readlines()]
        limits = list()
        start_node  = list()
        #coordinate = np.matrix((len(row),3))
        X = list()
        Y = list()
        Z = list()
        start_blade_index = list()
        number_blade      = list()
        blade_id_storage  = list()

        """ PART 1
            storage of:
            1) first line to be read for the blade indexes
            2) number of the blade indexes
        """
        for j in range(len(row)):

            if 'MARKER_TAG' and surface in row[j]:
                start_blade_index.append(j)                   # start read the blade index --> 
                integer = int(filter(str.isdigit, row[j+1]))  # how many index are there into the blade mesh?
                number_blade.append(integer)                  # storage of the integer into an empty list
    
        start = np.array(start_blade_index)
        start = int(start+1 )
        end   = np.array(number_blade)
        end   = int(end)
        end   = int(end + start +1)

        """ PART 2
            storage of the blade indexes:
        """
        for j in range(len(row)):
            if (j > start) and (j < end ):
                elem = row[j].split()
                blade_id_storage.append(elem[1])

        """ PART 3  
            check of mesh nodes (total number!)
        """
        for j in range(len(row)):

            if 'NPOIN' in row[j]:
                upper_limit = int(j) 
                limits.append(upper_limit)
                start_value = int(len(limits)-1)
                start_node.append(start_value)
                # added on the 19/12/2018: the first column should be selected.
                number_of_points = row[j].split()
                points = number_of_points[1]
                #print 'from line 71 of the for statement: the number of points should be:', points
                #
                #number_points =  int(filter(str.isdigit, row[j]))
                limits.append(points)
                #print 'the number of nodes inside the mesh is:', number_points
                index_of_last_pnt = int(limits[0]) + int(limits[1])
                limits.append(index_of_last_pnt)
 
        #print 'the original list is:',start_node
        indx  =  (np.array(start_node))
        indx  =  int(indx) # INDEX OF LIMITS IN WHICH THE FIRST LINE TO BE READ IS STORED.
        #print 'cast operation over index:', indx
        start = limits[indx]

        print 'the mesh will be read for indexes bigger than (vim index)', start+1
        print 'limits:', limits
        print 'the last line (index of vim ) will be smaller than:', (limits[2]+2)
        for j in range(len(row)):    
            if (j > start) and (j < (limits[2]+1)):
               elem = row[j].split()
               X.append( elem[0])
               Y.append( elem[1])
               last = int(len(elem) - 1)
               Z.append( elem[last]) # here we are in two dimensions! so it will be be the INDEX array !

        x = np.array([X]).T
        y = np.array([Y]).T
        z = np.array([Z]).T

        #print 'the final x is:', x                         #ok
        coord = np.matrix(np.hstack((z,x,y)))
        #print 'after stack:', coord                        #ok
        #print 'the size of coordinates is:', (coord.shape) #ok

        """ PART 4 :
            comparison between index of blade and index of total nodes
        """ 
        sort_blade_X  = list()
        sort_blade_Y  = list()

        for i in range(len(coord)):
            for j in range(len(blade_id_storage)):
                if coord[i,0] == blade_id_storage[j]:
                    sort_blade_X.append(coord[i,1])
                    sort_blade_Y.append(coord[i,2])

        #print 'the sorted X are:', sort_blade_X
        #print 'is the length of X equal to 192?', len(sort_blade_X) # ok !

        """ PART 5:
            storage of the blade coordinates inside an external file
        """
        xx = np.array([sort_blade_X]).T
        yy = np.array([sort_blade_Y]).T
        ii = np.array([blade_id_storage]).T

        #print 'the final x is ', xx
        #matrix = np.hstack((xx,yy,ii))
        matrix = np.hstack((ii,xx,yy))

        with open('sorted_blade.txt', 'w') as blade:
            blade.write(str(matrix))

        with open('sorted_blade.txt', 'r') as r:
            with open('blade_points.txt', 'w') as w:
                data = r.read()
                data = data.replace("'", "")
                data = data.replace("]", "")
                data = data.replace("[", "")
                w.write(data)
        os.system('rm sorted_blade.txt')
        print '--------------------------------------------------------------'
        return matrix

    def TransformMesh(self, surface, scale=None, translate=None, rotate=None):
        """ This method has been written to transform the mesh:
            the inputs for the manipulation of the geometry are optional, 
            while the parameter 'surface' is compulsory.

            a scale action is applied with a constant factor 'scale'
            a translation action is imposed through an array [x,y,z]
            a rotation matrix is imposed by specifing the angle in degree.
        """
        mesh = self.ExtractSurface(surface)
        coordinates = np.matrix(mesh)
        #coordinates = coordinates.astype(float)
        #new_mesh = np.zeros(mesh.shape)
        new_mesh_0 = np.zeros(len(mesh))
        new_mesh_1 = np.zeros(len(mesh))
        new_mesh_2 = list()
        
        if scale is None:
            scale_factor = 1.
        else:
            scale_factor = scale
        if translate is None:
            displacement = np.array([0., 0.])
        else:
            displacement = np.array(translate)
        if rotate is None:
            angle = 0.
        else:
            angle = math.radians(rotate)
        
        # 1) translate the old mesh
        # 2) scale factor applied to all the coordinates
        for i in range(len(coordinates)):
            #new_mesh[i,1] = (float(coordinates[i,1]) + float(displacement[0])) * scale_factor
            #new_mesh[i,2] = (float(coordinates[i,2]) + float(displacement[1])) * scale_factor
            #new_mesh[i,0] = (int(coordinates[i,0]))
            #print 'added:',int(new_mesh[i,0])
            new_mesh_0[i] = (float(coordinates[i,1]) + float(displacement[0])) * scale_factor
            new_mesh_1[i] = (float(coordinates[i,2]) + float(displacement[1])) * scale_factor
            new_mesh_2.append( str(int(coordinates[i,0])))

            # 3) rotation around the Z axis which direction is along the monitor entering vector
            #new_mesh[i,1] = new_mesh[i,1]*math.cos(angle) - new_mesh[i,2]*math.sin(angle)
            #new_mesh[i,2] = new_mesh[i,1]*math.sin(angle) + new_mesh[i,2]*math.cos(angle)

            new_mesh_0[i] = new_mesh_0[i]*math.cos(angle) - new_mesh_0[i]*math.sin(angle)
            new_mesh_1[i] = new_mesh_1[i]*math.sin(angle) + new_mesh_1[i]*math.cos(angle)
        new_mesh_0 = np.array([new_mesh_0]).T
        new_mesh_1 = np.array([new_mesh_1]).T
        new_mesh_2 = np.array([new_mesh_2]).T
        #new_mesh_2 = new_mesh_2.T
        new_mesh = np.hstack((new_mesh_2, new_mesh_0, new_mesh_1))
        print 'the new matrix is:', new_mesh

        with open('transf_mesh.txt', 'w') as draft:
            draft.write(str(new_mesh))
        with open('transf_mesh.txt', 'r') as draft:
            with open('transformed_mesh.txt', 'w') as blade:
                data = draft.read()
                data = data.replace("]", "")
                data = data.replace("[", "")
                data = data.replace("'", "")
                blade.write(data)
        os.system('rm transf_mesh.txt')
            
        return new_mesh

    def SortedBlade(self, surface):
        """ This method has been written to associate a curvilinear coordinate to each point of the blade surface.
            The blade points which belong to the return statement of 'ExtractSurface' are sorted from the
            leading edge (with respect of the hypothesis that this section corresponds to the origin of 
            the axis).
        """
        blade = self.ExtractSurface(surface)
        sorted_blade = blade

        distance_index = list()
        distance_array = list()
        # 1 ) find the leading edge: the geometrical origin of the axis;
        x = np.array(blade[:,1])
        y = np.array(blade[:,2])
        x = (x.astype(np.float))
        y = (y.astype(np.float))

        minimum_x = float(min(x))
        minimum_y = float(min(y))
        #print 'the minimum values are:', minimum_x, minimum_y
        for i in range(len(blade)):
            actual_x = np.array(blade[i,1])
            actual_x = float(actual_x.astype(np.float))
            actual_y = np.array(blade[i,2])
            actual_y = float(actual_y.astype(np.float))
            #print 'type of element inside blade:', type(actual_x), 'the value is:', actual_x
            #print 'while the checking type is:', type(minimum_x), 'the value is:', minimum_x

            if (actual_x == minimum_x): 
                distance_index.append(i)
                #sorted_blade.append(blade[i,0:2])
            #if (actual_y == minimum_y ):
            #    print 'for th y array it has been found:', i
       
        print '\n------------------------------------remember that:---------------------------------------------------------'
        print '\nin the present loop the distinctive factor to find out the geometrical leading edge has been assigned to '
        print '\nthe minimum value of X. Y is a consequence;'
        print '\n(see for statement between lines 232-237 of ExtractCurvCoord method inside meshDeform class. )'
        #print '\nthe coordinates of the sorted blade are:', sorted_blade
        print '\n-----------------------------------------------------------------------------------------------------------'
        
        # 2 ) find the distances from the first point of the list          
        #print 'the distance_array is:', distance_index

        leading_edge_index = int(distance_index[0])
             
        x_ref = float(blade[leading_edge_index, 1])
        y_ref = float(blade[leading_edge_index, 2])


        middle = len(sorted_blade) +1
        sorted_blade1 = sorted_blade[leading_edge_index:middle]
        sorted_blade2 = sorted_blade[0:(leading_edge_index)]
        blade_from_le = np.vstack((sorted_blade1, sorted_blade2))

        ##print 'length of sorted_blade:', len(sorted_blade)
        with open('blade_from_le_0.txt', 'w') as b :
            b.write(str(blade_from_le))
        with open('blade_from_le_0.txt', 'r') as r:
            with open('blade_from_le.txt', 'w') as b:
                data = r.read()
                data = data.replace("[", "")
                data = data.replace("]", "")
                data = data.replace("'", "")
                data = data.replace("array(", "")
                data = data.replace(")", "")
                data = data.replace("dtype=|S21", "")
                data = data.replace(",", "")
                b.write(data)
        os.system("rm blade_from_le_0.txt")

        return blade_from_le
    
    def ApplyBump(self, bump_array=None, surface=None):
        """ The aim of this method is to modify the blade geometry:
            starting from a surface in which all the points are 
            sorted from the leading edge, the curvilinear coordinate
            is calculated, the results are non-dimensionalized and
            the Hicks-Henne operator is applied over this new
            reference system.
        """
        
        if surface is None:
            print 'A surface name must be given as input to the method ApplyBump'
        else:
            srf = surface
        if bump_array is None:
            print 'A Bump matrix must be given as input to the method ApplyBump'
        else:
            b = np.matrix(bump_array)
 
        air = np.matrix(self.SortedBlade(srf))
        with open('bump.txt', 'w') as bw:
            bw.write(str(b))
        with open('bump.txt', 'r') as r :
            with open('bump_new.txt', 'w') as new:
                data = r.read()
                data = data.replace("[", "")
                data = data.replace("]", "")
                new.write(str(data))
        os.system("rm bump.txt")
        os.system("cp DVPerturbation/newhickshenne.m .")
        os.system("./newhickshenne.m")
        os.system("rm newhickshenne.m")
        os.system("rm bump_new.txt")
       

        """
        # air contains the blade coordinates, sorted from the geometrical leading edge.
        air = np.matrix(self.SortedBlade(srf))
        n   = len(air)
        # into the variable sc it will be stored the number of the rows inside the bump matrix
        sc  = len(b)
        # t2 is one the parameters to determine the shape of the HicksHenne curve.
        t2  = 6.65
        
        columns = int(air.size/len(air))  # how many columns are there inside the air matrix?
        c = int(columns + (7 + sc))
        tmp = np.zeros((len(air),c))      # temporary results matrix

        for i in range(n-1):
            tmp[i,0] = int(air[i,0])
            tmp[i,1] = air[i,1]
            tmp[i,2] = air[i,2]

            ii = int(i +1)

            # the following line calculates the distance between neighbouring points along the surface
            tmp[ii,3] = np.sqrt((tmp[ii,1]-tmp[i,1])**2 + (tmp[ii,2]-tmp[i,2])**2)
            # cumulative curvilinear coordinate:
            tmp[ii,4] = tmp[i,4] + tmp[ii,3]
        # non-dimensional cumulative curvilinear coordianate:
        tmp[:,5] = tmp[:,4]/tmp[n-1,4]
        for j in range(sc):
            jmin = int(j - 1)
            jmax = int(j + 1)
            if (j == 0 ):
                jmin = int(sc-1)
            if (j == sc-1):
                jmax = 0
            scenter = b[j,0]
            smin    = b[jmin,0]
            smax    = b[jmax,0]
            amp     = b[j,1]
            diffmin = scenter - smin
            diffmax = smax - scenter
            if (diffmin < 0):
                diffmin = diffmin + 1
                smin    = smin - 1
            if (diffmax < 0):
                diffmax = diffmax + 1
                smax    = smax + 1
            for i in range(n):
                s = tmp[i,5]
                if (s - smin > 1.):
                    s = s - 1.
                if (smax - s > 1.):
                    s = s + 1.
                if (s >= smin) and (s < scenter):
                    # first half of HH bump function on control point j starting from zero point on point j-1
                    tmp[i,5+j] = amp*(np.sin(.5*np.pi*(s-smin)/diffmin))**t2
                if (s >= scenter) and (s<= smax):
                    # second half of the HH bump function on the control point j reaching zero on point j+1
                    tmp[i,5+j] = amp*(np.sin(.5*np.pi*(smax-s)/diffmax))**t2
                else:
                    tmp[i,5] = 0.

        for j in range(sc):
            for i in range(n):
                im = int(i-1)
                ip = int(i+1)
                if (i==0):
                    im = int(n-1)
                if (i== n-1):
                    ip = 0
                xb = tmp[im,1]
                xm = tmp[i,1]
                xf = tmp[ip,1]
                #
                yb = tmp[im,2]
                ym = tmp[i,2]     
                yf = tmp[ip,2]
                #
                Dy = yf - yb
                Dx = xf - xb
                #
                M = -Dx/Dy
                a = np.arctan(M)
                d = tmp[i,5+j]
                if (Dx>0 and Dy<0):
                    # first quarter
                    xsign = np.sign(b[j,1])
                    ysign = xsign
                if (Dx<0) and (Dy<0):
                    # second quarter
                    xsign = np.sign(b[j,1])
                    ysign = -xsign
                if (Dx<0) and (Dy>0):
                    # third quarter
                    xsign = -np.sign(b[j,1])
                    ysign = xsign
                if (Dx>0) and (Dy>0):
                    # fourth quarter
                    xsign = -np.sign(b[j,1])
                    ysign = -xsign
                #print 'to verify: from line 408 of apply_bump method, xsign and ysign should be stored', xsign, ysign # ok --> print values to output
                tmp[i,6+sc] = tmp[i,6+sc]+ xsign*(abs(d*np.cos(np.arctan(M))))
                tmp[i,7+sc] = tmp[i,7+sc]+ ysign*(abs(d*np.sin(np.arctan(M))))
                #print 'from line 411 of apply_bump method: the ith element inside the matrix is', tmp[i,7+sc]
        tmp[:,8+sc]=tmp[:,6+sc] + tmp[:,1]
        tmp[:,9+sc]=tmp[:,7+sc] + tmp[:,2]
        
        # ----------------------------------------------------------------#
        # output section: which lines should be exported?
        # ----------------------------------------------------------------#
        #a  = (np.array([tmp[:,0]]).T)
        a = list()
        for i in range(len(tmp)):
            #a[i] = str(int(a[i]))
            a.append(str(int(tmp[i,0])))
            #print 'from for statement: the ith element of a is', a[i]
        a = np.array([a]).T
        #print 'len of a:', len(a)
        #print 'the size of a is:', a.shape
        b = np.array([tmp[:,8+sc]]).T
        #print 'len of b:', len(b)
        c = np.array([tmp[:,9+sc]]).T
        sp = np.hstack((a, b, c))
        #print sp
        #print 'the type of the final matris is :' , type(sp)

        with open('apply_bump_0.txt','w') as first:
            first.write(str(sp))
        with open('apply_bump_0.txt','r') as first:
            with open('appy_bump.txt','w') as sec:
                data = first.read()
                data = data.replace("[","")
                data = data.replace("]","")
                data = data.replace("'","")
                sec.write(data)
        os.system("rm apply_bump_0.txt")
        #-----------------------------------------------#
        with open('tmp.txt','w') as t:
            t.write(str(tmp))
        #-----------------------------------------------#
        #print tmp
        #print 'the lenght of tmp matrix is:', tmp.shape
        return tmp
        """
        #print 'end from applybump'

    def ApplyDEF(self, config_file=None):
        """ this method aims to apply the DEF module on the original mesh 
        """
        if config_file is None:
            print 'from ApplyDEF method int omeshDeform class: the name of the configuration file MUST BE GIVEN'
        else:
            os.system("SU2_DEF "+config_file)


    def VerifyIntegrity(self, mesh_new=None, mesh_old=None):
      """ This method has been written to check the maximum displacement of the nodes
          after the usage of the DEF module.
      """
      def nodes_filter(name,number):
        filename = name
        #filename = mesh_old
                                                                                                                     
        row = [x for x in open(filename).readlines()]
        limits = list()
        start_node  = list()
        #coordinate = np.matrix((len(row),3))
        X = list()
        Y = list()
        Z = list()
        start_blade_index = list()
        number_blade      = list()
        blade_id_storage  = list()
                                                                                                              
        """ PART 3  
            check of mesh nodes (total number!)
        """
        for j in range(len(row)):
                                                                                                              
            if 'NPOIN' in row[j]:
                upper_limit = int(j) 
                limits.append(upper_limit)
                start_value = int(len(limits)-1)
                start_node.append(start_value)
                # added on the 19/12/2018: the first column should be selected.
                number_of_points = row[j].split()
                points = number_of_points[1]
                #print 'from line 71 of the for statement: the number of points should be:', points
                #
                #number_points =  int(filter(str.isdigit, row[j]))
                limits.append(points)
                #print 'the number of nodes inside the mesh is:', number_points
                index_of_last_pnt = int(limits[0]) + int(limits[1])
                limits.append(index_of_last_pnt)
                                                                                                              
        #print 'the original list is:',start_node
        indx  =  (np.array(start_node))
        indx  =  int(indx) # INDEX OF LIMITS IN WHICH THE FIRST LINE TO BE READ IS STORED.
        #print 'cast operation over index:', indx
        start = limits[indx]
                                                                                                              
        print 'the mesh will be read for indexes bigger than (vim index)', start+1
        print 'limits:', limits
        print 'the last line (index of vim ) will be smaller than:', (limits[2]+2)
        for j in range(len(row)):    
            if (j > start) and (j < (limits[2]+1)):
               elem = row[j].split()
               X.append( float(elem[0]))
               Y.append( float(elem[1]))
               last = int(len(elem) - 1)
               Z.append( elem[last]) # here we are in two dimensions! so it will be be the INDEX array !
        
                                                                                                              
        x = np.array([X]).T
        y = np.array([Y]).T
        z = np.array([Z]).T
                                                                                                              
        #print 'the final x is:', x                         #ok
        coord = np.matrix(np.hstack((z,x,y)))
        np.set_printoptions(threshold=np.nan)

        with open('draft_verify.txt', 'w') as draft:
                draft.write(str(coord))
        with open('draft_verify.txt', 'r') as draft:
            names = str('mesh_verify'+str(number)+'.txt')
            with open(names, 'w') as mesh:
                data = draft.read()
                data = data.replace("[", "")
                data = data.replace("]", "")
                data = data.replace("'", "")
                mesh.write(str(data))
        os.system("rm draft_verify.txt")
        
      #------------------------------------------------------------
      # part 2: reading of the new mesh after the ApplyDEF method
      #------------------------------------------------------------
      array_of_mesh = np.array([mesh_old, mesh_new])
      #print 'the length of arrau is:', len(array_of_mesh)
      for j in range(len(array_of_mesh)):
          #print 'from for statement: the value of j is', j
          nodes_filter(array_of_mesh[j],j)
      #------------------------------------------------------------
      # part 3: comparison between meshes
      #------------------------------------------------------------




      #return coord

