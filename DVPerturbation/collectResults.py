""" This class has been written to automatically collect all the methods results inside a subdirectory,
    together with a readme fle, which contains the list of the methods involved.
"""
import numpy as np
import os
import datetime

""" This class takes as inputs:
    
    : class_name       = string, it is the name of the class
    : involved_methods = string, it contains the names of all the methods that have been invoked

"""

class CollectResults(object) :
    """ ATTRIBUTES
    """
    def __init__(self, class_name=None, involved_outputs=None, invoked_method=None, description=None):
        print 'Collecting outputs inside sub-directories'
        self.class_name = class_name
        self.involved_outputs = involved_outputs
        self.invoked_method   = invoked_method
        self.description = description

    """ METHODS
    """
    def Collect(self):
        """ This Method collect the results of the run of a method into a dedicated sub-directory
        """
        new_dir = str(self.class_name+"_OUTPUTs")
        sub_new = str(self.invoked_method+"_OUTPUTs")
        os.system("mkdir "+sub_new)
        
        for i in range(len(self.involved_outputs)):
            os.system("cp "+ self.involved_outputs[i] +" "+sub_new+"/")
        readme_content = str('This subdirectory contains the output of the class '+ self.class_name+'\n The method used is: '+self.invoked_method+'\n Description of the outputs: '+self.description)
        with open('LOGFILE_README.txt', 'w') as note: 
            note.write(str(readme_content))
        os.system("mv LOGFILE_README.txt "+sub_new+"/")
        os.system("mv "+sub_new+"/ "+new_dir+"/")

        """ In the following lines is reported the list of command that should be given 
            inside each class to generate the dedicated subfolder.
            to add:

            1) at the beginning of the class:
            from collectResults import CollectResults
            
            2) inside the __init__ method:
            self.class_name = 'class_name'
            os.system("mkdir "+str(self.class_name+"_OUTPUTs"))
            
            3) at the end of each method:
            # part for saving the results inside a sub-directory.                
            description=' (...) '
            involved_outputs=np.array(['file1.txt', ...,'filen,txt'])
            # only if the output is recursive:
            #for i in range(5):
            #   file_ith = str('-r CONFIG'+str(i)) 
            #   involved_outputs.append(file_ith)
            #   print 'the ith name is:', involved_outputs[i]
            invoked_method = 'SortedBlade'
            collection = CollectResults(class_name=self.class_name, involved_outputs=involved_outputs, invoked_method=invoked_method, description=description)
            collection.Collect()
        """
    def CheckWritingDir(self):
        """ This method has been writen to control the directory before a directory would be created
        """
        now = datetime.datetime.now()
        hour = str(str(now.hour)+str(now.minute)+str(now.second))
        list_of_files =  [f for f in os.listdir('.')]
        if str(self.class_name+'_OUTPUTs') in list_of_files:
            os.system('mv '+str(self.class_name+'_OUTPUTs')+' '+self.class_name+'_OUTPUTs'+str(hour))
        #print 'from the collect class: the list of files is:', list_of_files
        #print 'and its type is:', type(list_of_files)



