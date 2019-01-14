""" This class has been written to automatically run the post processing phase
"""
import numpy as np
import os

""" This class takes as inputs:
    
    : file_name = string, it specifies the file that will be processed
    : aim       = string, it determines the specific aim of the calculations
                  for example: 'gradients'
"""
class PostPro(object):
    """ ATTRIBUTES
    """
    def __init__(self):
        print 'Beginning of the Post Processing phase ..'

        
    """ METHODS
    """
    def Gradient(self, file_name=None):
        """ This method calculates the gradient starting from a file:
            this later must be into an authorized list 
        """
        authorized  = ['surface_sens.dat']
        dimension   = list() # this element contains the dimension of the problem
        net_results = list() # in this list it will be stored the indexes of the net matrix that will be processed

        if file_name in authorized:
            #
            data    = [x for x in open(file_name).readlines()]
            #results = np.zeros((len(data),2))     # final matrix of results
            if (file_name == 'surface_sens.dat'):
                    for i in range(len(data)):
                        row = data[i].split()
                        if 'VARIABLES' in data[i]:
                            if 'z' in row[i]:
                                dimension.append(3)
                            else:
                                dimension.append(2)
                        
                        if 'ZONE' in data[i]:
                            first_net = (i)
                            net_results.append(first_net)

                            check_number = row[2]
                            net_results.append(int(filter(str.isdigit ,check_number)))
                        
                        #print 'check of the net_results list inside a for statement maybe the list will be visible after its end:', net_results

                    results = np.zeros(( int(net_results[1]) , int(int(dimension[0]+1) )))     # final matrix of results
                    for i in range(len(data)):
                        row = data[i].split()
                        if (i > int(net_results[0])) and (i<= int(int(net_results[0]) + int(net_results[1]))):
                            results[int(i-(int(net_results[0]+1))),0] = row[0] # X
                            results[int(i-(int(net_results[0]+1))),1] = row[1] # Y
                            if dimension[0]==3:
                                results[int(i-int(net_results[0]+1)), 2] = row[2] # Z
                            results[int(i-(int(net_results[0]+1))),int(dimension[0])] = row[int(dimension[0]*2)] # sensitivities
                            #results[int(i-(int(net_results[0]+1))),int(dimension[0]+1)] = i                      # fictitious index
                        
                    # the matrix whixh contains the X, Y coordinates in defined
                    #-------------------------------------------------------------------------------------------------------------------#
                    # in the following lines the algorithm to sort the points inside the blade has been implemented
                    # because the same logic has been used to organize the TMP matrix of results: the corrispondence 
                    # between coordinates and sensitivities is possible through the sorting calculation
                    # the last columns contains the sensitivities.

                    blade = results

                    #print 'the initial blade file is', blade
                    distance_index = list()
                    distance_array = list()
                    # 1 ) find the leading edge: the geometrical origin of the axis;
                    x = np.array(blade[:,0]).astype(np.float)
                    y = np.array(blade[:,1]).astype(np.float)

                    minimum_x = float(min(x))
                    minimum_y = float(min(y))
                    #print 'the minimum values are:', minimum_x, minimum_y
                    for i in range(len(blade)):
                        actual_x = np.array(blade[i,0])
                        actual_x = float(actual_x.astype(np.float))
                        actual_y = np.array(blade[i,1])
                        actual_y = float(actual_y.astype(np.float))

                        if (actual_x == minimum_x):
                            distance_index.append(i)
                            #sorted_blade.append(blade[i,0:2])
                        #if (actual_y == minimum_y ):
                        #    print 'for th y array it has been found:', i

                    print '\n------------------------------------remember that:---------------------------------------------------------'
                    print '\nThe method Gradient inside PostPro class is working:'
                    print '\nin the present loop the distinctive factor to find out the geometrical leading edge has been assigned to '
                    print '\nthe minimum value of X. Y is a consequence;'
                    #print '\nthe coordinates of the sorted blade are:', sorted_blade
                    print '\n-----------------------------------------------------------------------------------------------------------'

                    # 2 ) find the distances from the first point of the list          
                    #print 'the distance_array is:', distance_index

                    leading_edge_index = int(distance_index[0])

                    middle = len(blade) +1
                    sorted_blade1 = blade[leading_edge_index:middle, :]
                    sorted_blade2 = blade[0:(leading_edge_index), :]
                    blade_from_le = np.vstack((sorted_blade1, sorted_blade2))                      
                    #print 'the sorted blade is :', blade_from_le
        #           #
        #----------------------------------------------------------------------------------------------------#
        else: 
            print 'the file name specified is not in the authorized list for the calculation of the gradient.'
            print 'see line 25 of the method Gradient inside the class postPro.'
            
        # part of gradient calculation starting from the sensitivity data
        # CONSTRUCTION OF THE MATRIX A
        
        rows_from_tmp   = [x for x in open('TMP.dat').readlines()]
        first_counter   = int(3+dimension[0])
        bump_numerosity = int(len(rows_from_tmp[0].split()) - (2*dimension[0] + dimension[0]+4) )
        bump_from_tmp = np.zeros((len(rows_from_tmp), bump_numerosity))
        #print 'from line 123 of postPro class: the bump numerosity is:', bump_numerosity
        for i in range(len(rows_from_tmp)):
            elem_from_tmp = rows_from_tmp[i].split()
            for j in range(bump_numerosity):
                # extract the bump matrix 
                index = int(dimension[0]+4 + j)
                bump_from_tmp[i,j] = elem_from_tmp[index]
        #print 'the bump matrix is:',  bump_from_tmp
        bump_from_tmp = bump_from_tmp.T
        #print 'it should be the transpose:', bump_from_tmp
        #print 'and the shape is:', bump_from_tmp.shape
        max_array = list() # it contains the maximum values of each row inside the bump_from_tmp matrix
        # original ------------------------------------------
        #for i in range(len(bump_from_tmp)):
        #    max_array.append(max(bump_from_tmp[i,:]))

        #max_array = np.array(max_array)
        # ---------------------------------------------------
        bump_entity = [x for x in open('bump_new.txt').readlines() ]
        #print 'this is the array with all the maximum values:', max_array
        non_dim_bump = bump_from_tmp
        columns = len(rows_from_tmp)
        for i in range(len(bump_from_tmp)):
            b_entity = bump_entity[i].split()
            for j in range(columns):
                #den = max_array[i]
                den = float(b_entity[1])
                #print 'this is inside  the for statement:line 150, it should be the bump amplitude',den
                #print 'the type of the denominator is:', type(den)
                if (den != 0.0) :
                    non_dim_bump[i,j] = non_dim_bump[i,j]/den
                else:
                    non_dim_bump[i,j] = 0.
        #print 'it should be a non dimensional (because the max =1) array:'
        #print non_dim_bump
        product = np.matmul(non_dim_bump, np.array([blade_from_le[:,int(dimension[0])]]).T)
        print 'the gradients are: -see lines 153 of Gradient method, inside PostPro class-', product
        #print 'the row of sensitivity is:', np.array([results[:,int(dimension[0])]]).T
        #print 'the shape of the column is:', np.array([results[:,int(dimension[0])]]).T.shape

        

        






