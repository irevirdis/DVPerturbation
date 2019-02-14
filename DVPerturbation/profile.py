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
class Profile(object):
    # attributes of the class
    
    def __init__(self, points=None, indices=None):
        self.class_name = 'Profile'
        self.points = points
        if (indices != None):
            self.indices=indices
        else:
            self.indices=np.array([[x for x in range(len(points))]]).T
            
        #os.system("mkdir "+str(self.class_name+"_OUTPUTs"))
    
    # methods
    
    def SortProfile(self):
        """ This method has been written to associate a curvilinear coordinate to each point of the blade surface.
            The blade points which belong to the return statement of 'ExtractSurface' are sorted from the
            leading edge (with respect of the hypothesis that this section corresponds to the origin of 
            the axis).
        """
        distance_index = list()
        distance_array = list()
        # 1 ) find the leading edge: the geometrical origin of the axis;
        x = np.array(self.points[:,0])
        y = np.array(self.points[:,1])

        x = (x.astype(np.float))
        y = (y.astype(np.float))

        minimum_x = float(min(x))
        minimum_y = float(min(y))
        #print 'the minimum values are:', minimum_x, minimum_y
        for i in range(len(self.points)):
            actual_x = np.array(self.points[i,0])
            actual_x = float(actual_x.astype(np.float))
            actual_y = np.array(self.points[i,1])
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
             
        x_ref = float(self.points[leading_edge_index, 0])
        y_ref = float(self.points[leading_edge_index, 1])


        middle = len(self.points) +1
        sorted_blade1 = self.points[leading_edge_index:middle, :]
        sorted_blade2 = self.points[0:(leading_edge_index), :]
        indices1 = self.indices[leading_edge_index:middle]
        indices2 = self.indices[0:(leading_edge_index)]
        blade_from_le = np.vstack((sorted_blade1, sorted_blade2))
        indices_from_le = np.vstack((indices1, indices2))

        ##print 'length of sorted_blade:', len(sorted_blade)
        np.set_printoptions(threshold=np.nan)
        with open('blade_from_le_0.txt', 'w') as b :
            b.write(str(np.hstack((indices_from_le,blade_from_le))))
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

        print blade_from_le[1,1], blade_from_le[0,1]

        if float(blade_from_le[1,1]) < float(blade_from_le[0,1]):
            blade_from_le = blade_from_le[::-1]
            indices_from_le = indices_from_le[::-1]
 
        
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
        self.sorted_points=blade_from_le
        self.sorted_indices=indices_from_le

    def ApplyBump(self, bump_array=None):
        """ The aim of this method is to modify the blade geometry:
            starting from a surface in which all the points are 
            sorted from the leading edge, the curvilinear coordinate
            is calculated, the results are non-dimensionalized and
            the Hicks-Henne operator is applied over this new
            reference system.
        """
        
        if bump_array is None:
            print 'A Bump matrix must be given as input to the method ApplyBump'
        else:
            b = np.matrix(bump_array)
 
        self.SortProfile()

        #with open('bump.txt', 'w') as bw:
        #    bw.write(str(b))
        #with open('bump.txt', 'r') as r :
        #    with open('bump_new.txt', 'w') as new:
        #        data = r.read()
        #        data = data.replace("[", "")
        #        data = data.replace("]", "")
        #        new.write(str(data))
        #os.system("rm bump.txt")
        #os.system("cp DVPerturbation/newhickshenne.m .")
        #os.system("./newhickshenne.m")
        #os.system("rm newhickshenne.m")
        #os.system("rm bump_new.txt")
        
        # sorted_points contains the blade coordinates, sorted from the geometrical leading edge.
        n   = len(self.sorted_points)
        # into the variable sc it will be stored the number of the rows inside the bump matrix
        sc  = len(b)
        # t2 is one the parameters to determine the shape of the HicksHenne curve.
        t2  = 6.65
        
        columns = int((self.sorted_points.size)/len(self.sorted_points))  # how many columns are there inside the sorted_points matrix?
        c = int(columns + (8 + sc))
        tmp = np.zeros((len(self.sorted_points),c))      # temporary results matrix

        for i in range(n):
            tmp[i,0] = int(self.sorted_indices[i])
            tmp[i,1] = self.sorted_points[i,0]
            tmp[i,2] = self.sorted_points[i,1]

        for i in range(1,n):
            # the following line calculates the distance between neighbouring points along the surface
            tmp[i,3] = np.sqrt((tmp[i,1]-tmp[i-1,1])**2 + (tmp[i,2]-tmp[i-1,2])**2)
            # cumulative curvilinear coordinate:
            tmp[i,4] = tmp[i-1,4] + tmp[i,3]

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
                    tmp[i,6+j] = amp*(np.sin(.5*np.pi*(s-smin)/diffmin))**t2
                elif (s >= scenter) and (s<= smax):
                    # second half of the HH bump function on the control point j reaching zero on point j+1
                    tmp[i,6+j] = amp*(np.sin(.5*np.pi*(smax-s)/diffmax))**t2
                else:
                    tmp[i,6+j] = 0.

        #for i in range(1,n):
            #print tmp[i,0], tmp[i,1], tmp[i,2], tmp[i,3], tmp[i,4], tmp[i,5], tmp[i,6], tmp[i,7], tmp[i,8], tmp[i,9]
        #print 'Done printing curvilinear coordinates'


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
                d = tmp[i,6+j]
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
                #print i, tmp[i,6+sc], tmp[i,7+sc]
                #print 'from line 411 of apply_bump method: the ith element inside the matrix is', tmp[i,7+sc]

        tmp[:,8+sc]=tmp[:,6+sc] + tmp[:,1]
        tmp[:,9+sc]=tmp[:,7+sc] + tmp[:,2]
        
        #for i in range(1,n):
        #    print tmp[i,0], tmp[i,1], tmp[i,2], tmp[i,3], tmp[i,4], tmp[i,5], tmp[i,6], tmp[i,7], tmp[i,8], tmp[i,9], tmp[i,10], tmp[i,11], tmp[i,12], tmp[i,13]
        #print 'Done printing curvilinear coordinates'


        self.deformed_points=np.hstack((np.array([tmp[:,8+sc]]).T, np.array([tmp[:,9+sc]]).T))    
        #sp = np.hstack((self.sorted_indices,  np.array([tmp[:,8+sc]]).T, np.array([tmp[:,9+sc]]).T))   

        """
        with open('apply_bump_0.txt','w') as first:
            #first.write(str(sp))
            first.write(str(self.sorted_indices))
        with open('apply_bump_0.txt','r') as first:
            with open('apply_bump.txt','w') as sec:
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


    def ExportDeformed(self,filename=None):
        if (filename == None):
            filename='surface_positions.dat'

        with open('sorted_positions.dat','w') as first:
            for i in range(len(self.sorted_indices)):
                first.write("%d %f %f \n" % (self.sorted_indices[i], self.sorted_points[i,0], self.sorted_points[i,1]))

        with open(filename,'w') as first:
            for i in range(len(self.sorted_indices)):
                first.write("%d %f %f \n" % (self.sorted_indices[i], self.deformed_points[i,0], self.deformed_points[i,1]))
        
        description='The output of this method is represented by the set of files bump_new.txt (which contains the matrix of applied bumps) and the files TMP.dat (the complete set of coordinates, applied bump and non-dimensional curvilinear coordinate) and '+filename+' (which contains the new coordinates X,Y of the blade surface)'
        involved_outputs= np.array([filename])   #list()
        invoked_method = 'ApplyBump'
        collection = CollectResults(class_name=self.class_name, involved_outputs=involved_outputs, invoked_method=invoked_method, description=description)
        collection.Collect()
    

    def Clockwise(self ):
        """ sort the points clockwise; this method has been written for 2D case;
            the input index is optional, used when the coordinates are associated
            to an index into a mesh.
        """ 
        matrix = self.points
        if self.indices is not None:
            index = self.indices
        copy_of_index = list()
        distance  = list()
        minimum_x = np.min(matrix[:,0])
        ref_first = list()
        x = list()
        y = list()

        for i in range(len(matrix)):
            x.append(matrix[i,0])
            if x[i] == minimum_x:
                ref_first.append(i)
            y.append(matrix[i,1])
            if index is not None:
                copy_of_index.append(index[i])
        ref_first = int(ref_first[0])

        sort_x = list()
        sort_y = list()
        sort_i = list()

        # initialize
        sort_x.append(x[ref_first])
        sort_y.append(y[ref_first])
        if index is not None:
            sort_i.append(copy_of_index[ref_first])
        ref_x = x[ref_first]
        ref_y = y[ref_first]
        
        x.pop(ref_first)
        y.pop(ref_first)
        
        time = len(x)
        #print 'the first value of sort_x is:', sort_x, 'and y:', sort_y
        for i in range(time):
            distance = list()            
            #print 'the list of remaining x:', x
            #print '--------------------- y:', y
            #print 'from external loop: the value of time is', i
            for j in range(time):
                #print 'from nested loop: time is',  j

                distance.append(np.sqrt((x[j]-ref_x)**2 + (y[j] - ref_y)**2))
   
            local_minimum = np.min(distance)
            idx = distance.index(local_minimum)
            #print 'the distances are:', distance
            #print 'the index of the smallest distance is', idx
            #print 'the present length of x is:', len(x)
            if len(x) != 0:
                sort_x.append(x[idx])   
                sort_y.append(y[idx])
            
                # new loop:
                ref_x = x[idx]
                ref_y = y[idx]

                distance.pop(idx)
                #print 'the length of distance is:', len(distance)
                x.pop(idx)
                y.pop(idx)
                if index is not None:
                    sort_i.append(copy_of_index[idx])
                    copy_of_index.pop(idx)
                time = len(x)
                #print 'the new value of time is:', time
                if time == 0:
                    break
                #print 'end of the nested-----------------------------------------'
            if time == 0:
                break
            #print 'end of the external ------------------------------------------'

        # clockwise:
        if (sort_x[0] < sort_x[1]):
            sort_x = sort_x[::-1]
            sort_y = sort_y[::-1]
            sort_i = sort_i[::-1]

        #print 'the points are:', sort_x
        with open('prova.txt', 'w') as w :                  
            for i in range(len(sort_x)):
                if sort_i is not None:
                    w.write(" %f  %f %d \n" % (sort_x[i], sort_y[i], sort_i[i]))               
                else:
                    w.write("%f %f \n" % (sort_x[i], sort_y[i]))

            
        



