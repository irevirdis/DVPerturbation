""" The aim of the present class is to deform the mesh nodes, starting from an su2 format.
"""
import numpy as np
import os
import math
from collectResults import CollectResults
from decimal import Decimal

""" This class takes as inputs:
    : param mesh_name = it is a string that will points the mesh file.

    : RETURNS: a new mesh file in which the nodes coordinates have been changed.
"""
class MeshDeform(object):
    # attributes of the class
    
    def __init__(self, mesh_name):
        self.mesh_name = mesh_name

        self.class_name = 'MeshDeform'
        obj = CollectResults(class_name=self.class_name)
        obj.CheckWritingDir()        
        self.indices=None
        self.points=None
        os.system("mkdir "+str(self.class_name+"_OUTPUTs"))
    
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
	I = list()
        start_blade_index = list()
        number_blade      = list()
        blade_id_storage  = list()

        """ PART 1
            storage of:
            1) first line to be read for the blade indexes
            2) number of the blade indexes
        """
        for j in range(len(row)):

            if 'MARKER_TAG' and surface in row[j].split():
                start_blade_index.append(j)                   # start read the blade index --> 
                integer = int(filter(str.isdigit, row[j+1]))  # how many index are there into the blade mesh?
                number_blade.append(integer)                  # storage of the integer into an empty list
                ## added on the 11th of Feb: check the presence of a word inside another one-----------#
                ## the following lines are an heuristics to check the structure of the row:
                print 'from line 63 of meshDeform: I found the wors inside the row:', row[j]
                elem = row[j].split()
                for t in range(len(elem)):
                    print 'from the heuristic: the element ',t,' is:', elem[t]
                #--------------------------------------HEURISTIC ---------------------------------------#
    
        start = np.array(start_blade_index)
        start = int(start+1 )
        end   = np.array(number_blade)
        end   = int(end)
        end   = int(end + start +1)

        """ PART 2
            storage of the blade indexes:
        """
        for j in range(len(row)):
            problem_dimension = len(row[j].split())
            #print 'this is the dimension of the problem: from line 68 of class meshDeform',problem_dimension
            if (j > start) and (j < end ) :
                if (problem_dimension==3):
                    elem = row[j].split()
                    blade_id_storage.append(elem[1])
                else:
                    elem = row[j].split()
                    for k in range(int(len(elem)-1)):
                    #print 'from line 73, the dimension is greater than 2D!'
                    # filter part: into blade_id_storage we do the stack of the indexes without repetions    
                        if elem[int(k+1)] not in blade_id_storage:
                            #print 'i am doing the storage of this element:', elem[int(k+1)]
                            blade_id_storage.append(elem[int(k+1)])
        np.set_printoptions(threshold=np.sys.maxsize)
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
        print 'from the line 108 of the class: the index array is :', indx, 'it should be a single value for each ZONE'
        #indx  =  int(indx) # INDEX OF LIMITS IN WHICH THE FIRST LINE TO BE READ IS STORED.
        #print 'cast operation over index:', indx
        # added on the 10th of february: adapt for NZONE different from 1:
        # indx has more than one element.
        for k in range(len(indx)):
            start = limits[int(indx[k])]

            print 'the mesh will be read for indexes bigger than (vim index)', start+1
            print 'limits:', limits
            print 'the last line (index of vim ) will be smaller than:', (limits[2]+2)
            for j in range(len(row)):    
                if (j > start) and (j < (limits[2]+1)):
                    elem = row[j].split()
                    X.append( elem[0])
                    Y.append( elem[1])
	            Z.append( elem[2])
                    last = int(len(elem) - 1)
                    I.append( elem[last]) # here we are in two dimensions! so it will be be the INDEX array !

            x = np.array([X]).T
            y = np.array([Y]).T
            z = np.array([Z]).T
	    i = np.array([I]).T

            #print 'the final x is:', x                         #ok
            coord = np.matrix(np.hstack((i,x,y,z)))
            #print 'after stack:', coord                        #ok
            #print 'the size of coordinates is:', (coord.shape) #ok

            """ PART 4 :
            comparison between index of blade and index of total nodes
            """ 
            sort_blade_X  = list()
            sort_blade_Y  = list()
            sort_blade_I  = list()
	    sort_blade_Z  = list()

            for j in range(len(blade_id_storage)):
                #print blade_id_storage[j]
                #print 'coords', coord[int(blade_id_storage[j]),0], blade_id_storage[j]
                nindex=int(blade_id_storage[j])#-1
                if (coord[nindex,0] == blade_id_storage[j]):
                    sort_blade_I.append(int(coord[nindex,0]))
                    sort_blade_X.append(round(Decimal(coord[nindex,1]),6))
                    sort_blade_Y.append(round(Decimal(coord[nindex,2]),6))
                    sort_blade_Z.append(round(Decimal(coord[nindex,3]),6))
                else:
                    print 'I am inside this'
                    for i in range(len(coord)):
                        if coord[i,0] == blade_id_storage[j]:
                            sort_blade_I.append(int(coord[i,0]))
                            sort_blade_X.append(round(Decimal(coord[i,1]),6))
                            sort_blade_Y.append(round(Decimal(coord[i,2]),6))
		            sort_blade_Z.append(round(Decimal(coord[i,3]),6))


            """ PART 5:
            storage of the blade coordinates inside an external file
            """

            xx = np.array([sort_blade_X]).T
            yy = np.array([sort_blade_Y]).T
	    zz = np.array([sort_blade_Z]).T
            ii = np.array([sort_blade_I]).T
            #--------------------------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # to be implemented only for 2D cases!
            #if float(yy[1,0]) < float(yy[0,0]):
            #    xx = xx[::-1]
            #    yy = yy[::-1]
            #    ii = ii[::-1]

            matrix = np.hstack((ii,xx,yy,zz))
            self.indices=np.array([sort_blade_I]).T
            self.points=np.hstack((xx, yy,zz))

        
            with open('sorted_blade_'+str(k)+'.txt', 'w') as blade:
                blade.write(str(matrix))

            with open('sorted_blade_'+str(k)+'.txt', 'r') as r:
                with open('blade_points_'+str(k)+'.txt', 'w') as w:
                    data = r.read()
                    data = data.replace("'", "")
                    data = data.replace("]", "")
                    data = data.replace("[", "")
                    w.write(data)
            os.system('rm sorted_blade_'+str(k)+'.txt')
            print '--------------------------------------------------------------'
        

        # part for saving the results inside a sub-directory.                
        description=' The output of the method is the file blade_points.txt which contains the coordinate of the blade, with the same order specified inside the original mesh file with SU2 format.'
        #involved_outputs='blade_points_'+k+'.txt'
        involved_outputs = list()
        for i in range(len(indx)):
            file_ith = str('blade_points_'+str(i)+'.txt') 
            involved_outputs.append(file_ith)
        #print 'the ith name is:', involved_outputs[i]
        invoked_method = 'ExtractSurface'
        collection = CollectResults(class_name=self.class_name, involved_outputs=involved_outputs, invoked_method=invoked_method, description=description)
        collection.Collect()
        global_res = list()
        global_res.append(matrix)
        #return matrix
        return global_res


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
        
        for i in range(len(coordinates)):
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
            
        

        # part for saving the results inside a sub-directory.                
        description=' The output of the method is a file transformed_mesh.txt, which contains only the coordiantes of the blade. after a rotation, translation and a scale action; if the parameter for this three actions are not given, the default calculation is rotation of 0 degree, translation of (0,0) and scale of 1.0; the angle of rotation should be expressed with degree; the translation is expressed with the coordinates (X,Y) of the displacement and the scale is a constant factor for resizing the blade.)'
        involved_outputs=np.array(['transormed_mesh.txt'])
        #for i in range(5):
        #file_ith = str('-r CONFIG'+str(i)) 
        #involved_outputs.append(file_ith)
        #print 'the ith name is:', involved_outputs[i]
        invoked_method = 'TransformMesh'
        collection = CollectResults(class_name=self.class_name, involved_outputs=involved_outputs, invoked_method=invoked_method, description=description)
        collection.Collect()
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
	z = np.array(blade[:,3])	

        x = (x.astype(np.float))
        y = (y.astype(np.float))
	z = (z.astype(np.float))

        #minimum_x = float(min(x))
        minimum_y = float(min(y))
	minimum_z = float(min(z))
        #print 'the minimum values are:', minimum_x, minimum_y
        for i in range(len(blade)):
            #actual_x = np.array(blade[i,1])
            #actual_x = float(actual_x.astype(np.float))
            actual_y = np.array(blade[i,2])
            actual_y = float(actual_y.astype(np.float))
            actual_z = np.array(blade[i,3])
            actual_z = float(actual_z.astype(np.float))
            #print 'type of element inside blade:', type(actual_x), 'the value is:', actual_x
            #print 'while the checking type is:', type(minimum_x), 'the value is:', minimum_x

            if (actual_z == minimum_z): 
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
	z_ref = float(blade[leading_edge_index, 3])


        middle = len(sorted_blade) +1
        sorted_blade1 = sorted_blade[leading_edge_index:middle, :]
        sorted_blade2 = sorted_blade[0:(leading_edge_index), :]
        blade_from_le = np.vstack((sorted_blade1, sorted_blade2))

        ##print 'length of sorted_blade:', len(sorted_blade)
        np.set_printoptions(threshold=np.sys.maxsize)
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

        
        # part for saving the results inside a sub-directory.                
        description=' The putput of the method is a file blade_from_le.txt, which contains only the coordiantes of the blade, sorted by the leading edge point (here the hypothesis is that the leading edge is on the point with minimum value of X)'
        involved_outputs=np.array(['blade_from_le.txt'])
        #for i in range(5):
        #file_ith = str('-r CONFIG'+str(i)) 
        #involved_outputs.append(file_ith)
        #print 'the ith name is:', involved_outputs[i]
        invoked_method = 'SortedBlade'
        collection = CollectResults(class_name=self.class_name, involved_outputs=involved_outputs, invoked_method=invoked_method, description=description)
        collection.Collect()
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
        #os.system("rm newhickshenne.m")
        #os.system("rm bump_new.txt")
       

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
        # part for saving the results inside a sub-directory.                
        description='The output of this method is represented by the set of files bump_new.txt (which contains the matrix of applied bumps) and the files TMP.dat (the complete set of coordinates, applied bump and non-dimensional curvilinear coordinate) and surface_positions.dat (which contains the new coordinates X,Y of the blade surface)'
        involved_outputs= np.array(['TMP.dat', 'surface_positions.dat'])   #list()
        #for i in range(2):
        #file_ith = str('-r CONFIG'+str(i)) 
        #involved_outputs.append(file_ith)
        #print 'the ith name is:', involved_outputs[i]
        invoked_method = 'ApplyBump'
        collection = CollectResults(class_name=self.class_name, involved_outputs=involved_outputs, invoked_method=invoked_method, description=description)
        collection.Collect()


    def ApplyDEF(self, config_file=None):
        """ this method aims to apply the DEF module on the original mesh 
        """
        if config_file is None:
            print 'from ApplyDEF method int omeshDeform class: the name of the configuration file MUST BE GIVEN'
        else:
            os.system("SU2_DEF "+config_file)

        # part for saving the results inside a sub-directory.                
        description='The output is the new mesh coordinates, after the application of the DEF module. The output is a .su2 mesh'
        involved_outputs= np.array(['mesh_out.su2']) #list()
        #for i in range(1):
        #file_ith = str('-r CONFIG'+str(i)) 
        #involved_outputs.append(file_ith)
        #print 'the ith name is:', involved_outputs[i]
        invoked_method = 'ApplyDEF'
        collection = CollectResults(class_name=self.class_name, involved_outputs=involved_outputs, invoked_method=invoked_method, description=description)
        collection.Collect()


    def VerifyIntegrity(self, mesh_new=None, mesh_old=None, bump_matrix=None):
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
                                                                                                              
        coord = np.matrix(np.hstack((z,x,y)))
        #np.set_printoptions(threshold=np.sys.maxsize)
        np.set_printoptions(thereshold=sys.maxsize)
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

      for j in range(len(array_of_mesh)):
          nodes_filter(array_of_mesh[j],j)
      #------------------------------------------------------------
      # part 3: comparison between meshes
      #------------------------------------------------------------
        
      # maximum bump entity
      b     = np.matrix(bump_matrix)
      max_b = 2*max(bump_matrix[:,1])

      row_old = [x for x in open('mesh_verify0.txt').readlines()]
      row_new = [y for y in open('mesh_verify1.txt').readlines()]
      
      mesh_corrected = np.zeros((len(row_new), 3))
      list_integer = list()
      
      for i in range(len(row_old)):
          old = row_old[i].split()
          new = row_new[i].split()
          difference_on_X = float(old[1]) - float(new[1])
          difference_on_Y = float(old[2]) - float(new[2])
          displacement = np.sqrt(difference_on_X**2 + difference_on_Y**2)
          if (displacement > abs(max_b)):
              list_integer.append(int(old[0]))
              
              mesh_corrected[i,0] = int(old[0])
              mesh_corrected[i,1] = float(old[1])
              mesh_corrected[i,2] = float(old[2])
          else:
              list_integer.append(str(int(new[0])))
              mesh_corrected[i,0] = int(new[0])
              mesh_corrected[i,1] = float(new[1])
              mesh_corrected[i,2] = float(new[2])

      a = np.array([list_integer]).T
      b = np.array([mesh_corrected[:,1]]).T
      c = np.array([mesh_corrected[:,2]]).T

      to_export = np.hstack((b,c,a))
      #----------------------------------------------------------------------------------
      #in the following lines only the coordinates of the corrected mesh will exported
      with open('correct2.txt', 'w') as corr:
            corr.write(str(to_export))
      with open('correct2.txt', 'r') as corr:
        with open('corrected_mesh.txt', 'w') as cor2:
            data = corr.read()
            #print 'from line 605 : the type of the data which is reading is:', type(data)
            data = data.replace("[","")
            data = data.replace("]","")
            data = data.replace("'","")
            cor2.write(data)
      os.system("rm correct2.txt")
      
      #----------------------------------------------------------------------------------
      # in this part a new mesh, with su2 format will be created.

      """ the aim of this section is to prepare a mesh file with a su2 format:
          the input 'old_mesh' will provide two main bloks inside the final 
          mesh; the input 'correct' will provide the block between NPOIN marker
          and NELEM marker.
      """
      mesh_old = self.mesh_name 
      mesh_rows = [x for x in open(mesh_old).readlines()]
      insert_correct = [y for y in open('corrected_mesh.txt').readlines()]

      central_part = list()

      for i in range(len(mesh_rows)):

          if 'NPOIN' in mesh_rows[i]:
              central_part.append(i)
              check_columns = mesh_rows[i].split()
              central_part.append(int(filter(str.isdigit, mesh_rows[i])))

      first_of_central =  float(central_part[0])
      last_of_central  =  float(central_part[1])

      first_block  = []

      for i in range(len(mesh_rows)):
          """
          this for statement end in the line containin the NPOIN marker:
          for this reason the i index is set on (first_central +1)
          """
          if (i <= first_of_central):
              with open('draft_to_full_mesh.txt', 'a') as draft:
                  draft.write(mesh_rows[i])
          if (i > int(first_of_central+1)) and (i <= int(last_of_central+2)):
              with open('draft_to_full_mesh.txt', 'a') as draft:
                  index = int(i-3 -last_of_central)
                  #--------------------------------------------------------------
                  # 11 / 01/ 2019
                  insert = insert_correct[index]
                  #print 'the type of my insert is:', type(insert)
                  insert_new = insert.lstrip()
                  draft.write(insert_new)
                  #--------------------------------------------------------------
                  #draft.write(insert_correct[index])

          if (i==int(last_of_central+2)):
              with open('draft_to_full_mesh.txt', 'a') as draft:
                  draft.write(str('\n'))
          if (i > int(last_of_central+1)):
              with open('draft_to_full_mesh.txt', 'a') as draft:
                  insert = mesh_rows[i]
                  insert_new = insert.lstrip()
                  draft.write(insert_new)

      #print 'the first of the central is :', first_of_central
      #print 'the last of the central is:', last_of_central

      os.system("mv draft_to_full_mesh.txt full_corrected_mesh.txt")
      os.system("dos2unix full_corrected_mesh.txt")
      os.system("cp full_corrected_mesh.txt mesh_verified.su2")

      # part for saving the results inside a sub-directory.                
      description='The output have been writen to check the maximum displacement of the nodes inside a mesh after the application of the DEF module.\n The file mesh_verify0.txt contains only the coordinates of the specified surface of the original blade; the file mesh_verify1.txt contains the coordinates of the blade surface after the application of the DEF module;\n the file corrected_mesh and full_corrected mesh.txt contain the new coordinates of  all the nodes (not only the blade) after the check of the maximum displacement (2*bump_entity); mesh_verified.su2 can be used for a simulation.'
      involved_outputs=np.array(['mesh_verify0.txt', 'mesh_verify1.txt', 'full_corrected_mesh.txt', 'mesh_verified.su2', 'corrected_mesh.txt'])
      #for i in range(5):
      #file_ith = str('-r CONFIG'+str(i)) 
      #involved_outputs.append(file_ith)
      #print 'the ith name is:', involved_outputs[i]
      invoked_method = 'VerifyIntegrity'
      collection = CollectResults(class_name=self.class_name, involved_outputs=involved_outputs, invoked_method=invoked_method, description=description)
      collection.Collect()


    def ExtractMesh(self):
        """ Method for the external export of the entire mesh (set of points from NPOIN zone)
        """
        dimension = list()
        index = list()
        mesh_file = [x for x in open(self.mesh_name).readlines()]
        for i in range(len(mesh_file)):
            if 'NPOIN' in mesh_file[i]:
                dimension.append(int(len(mesh_file[i+1].split()) -1 ))
                index.append(i)
                index.append(int(filter(str.isdigit, mesh_file[i])))
        nzone = int(len(index)/2)
        for j in range(nzone):
            read_the_index_1 = int(j*2)
            read_the_index_2 = int(read_the_index_1 +1)
            if dimension[0] == 2:
                matrix = np.zeros((int(index[read_the_index_2]), 2))
            else:
                matrix = np.zeros((int(index[read_the_index_2]), 3))

            for i in range(len(mesh_file)):
                if (i>index[read_the_index_1]) and (i<=(index[read_the_index_2] + index[read_the_index_1])):
                    position = int(i - index[read_the_index_1] -1)
                    elem = mesh_file[i].split()
                    matrix[position,0] = elem[0]
                    matrix[position,1] = elem[1]
                    if dimension[0] == 3 :
                        matrix[position,2] = elem[2]
            file_name = str('mesh_zone_'+str(j)+'.txt')
            with open(file_name, 'w') as w :
                for i in range(len(matrix)):
                    if dimension[0] == 2:
                        w.write("%f %f \n" % (matrix[i,0], matrix[i,1]))
                    else:
                        w.write("%f %f \n" % (matrix[i,0], matrix[i,1], matrix[i,2]))
        description = 'the output of this method is represented by the set of the points inside each NPOIN zone of the original mesh'
        involved_outputs = list()
        for i in range(nzone):
            to_save = str('mesh_zone_'+str(i)+'.txt')
            involved_outputs.append(to_save)
        invoked_method = 'ExtractMesh'
        collection = CollectResults(class_name = self.class_name, involved_outputs = involved_outputs, invoked_method = invoked_method, description= description)
        collection.Collect()


