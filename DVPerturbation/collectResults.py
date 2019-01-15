""" This class has been written to automatically collect all the methods results inside a subdirectory,
    together with a readme fle, which contains the list of the methods involved.
"""
import numpy as np
import os

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


