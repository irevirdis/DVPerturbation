""" This class represemts a writing generator of configuration files for the SU2 simulations.
"""
import numpy as np
import os
from collectResults import CollectResults

"""     this class takes as inputs:
	: param bump   =   it an array which contains the bump entity (DV variable inside a config file)
        : param pos    =   is an array which contains the position of each bump 
	:RETURNS: configuration_file.cfg con numerosita' pari a simulations.
"""
class CfgGenerator(object):
	# attributes will be specified in the future
	def __init__(self, bump=None, pos=None):
                print '---------------------------------------'
                print 'starting calculation from cfgClass:'
		self.n = len(bump)
                print 'the number of configurations is', self.n
		self.B = bump
		print 'the input B is:', self.B
		self.p = pos
		print 'the positions are', self.p
		#print 'from the init of cgf generator, the number of simulations will be:', self.n
                print '----------------------------------------'
                self.class_name = 'CfgGenerator'
                os.system("mkdir "+str(self.class_name+"_OUTPUTs"))

	#-------------------------------------------#
	# METHODS
	#-------------------------------------------#
	def ScanFolder(self):
		""" This method specifies the working directory in which configuration files will be written
		"""
		pwd = os.getcwd() 
		#print pwd
                files = os.listdir(pwd)
                return files, pwd

	def WriteDraft(self):
		""" this method will write inside the actual directory 
		"""		
		config='PHYSICAL_PROBLEM= NAVIER_STOKES\nKIND_TURB_MODEL= SA\nMATH_PROBLEM= DIRECT\nRESTART_SOL= NO\nMACH_NUMBER= 0.729\nAOA= 2.31\nFREESTREAM_TEMPERATURE= 288.15\nREYNOLDS_NUMBER= 6.5E6\nREYNOLDS_LENGTH= 1.0\nFIXED_CL_MODE= NO\nTARGET_CL= 0.724\nDCL_DALPHA= 0.2\nUPDATE_ALPHA= 5\nITER_DCL_DALPHA= 500\nEVAL_DOF_DCX= NO\nREF_ORIGIN_MOMENT_X= 0.25\nREF_ORIGIN_MOMENT_Y= 0.00\nREF_ORIGIN_MOMENT_Z= 0.00\nREF_LENGTH= 1.0\nREF_AREA= 1.0\nMARKER_HEATFLUX= ( AIRFOIL, 0.0 )\nMARKER_FAR= ( FARFIELD )\nMARKER_PLOTTING= ( AIRFOIL )\nMARKER_MONITORING= ( AIRFOIL )\nNUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES\nCFL_NUMBER= 10.0\nEXT_ITER= 50\nOBJECTIVE_FUNCTION= DRAG\nMUSCL_FLOW= YES\nSLOPE_LIMITER_FLOW= VENKATAKRISHNAN\nMUSCL_TURB= NO\nSLOPE_LIMITER_TURB= VENKATAKRISHNAN\nMUSCL_ADJFLOW= YES\nSLOPE_LIMITER_ADJFLOW= VENKATAKRISHNAN\nMUSCL_ADJTURB= NO\nSLOPE_LIMITER_ADJTURB= VENKATAKRISHNAN\nVENKAT_LIMITER_COEFF= 0.05\nADJ_SHARP_LIMITER_COEFF= 3.0\nLIMITER_ITER= 999999\nLAX_SENSOR_COEFF= 0.15\nJST_SENSOR_COEFF= ( 0.5, 0.02 )\nADJ_LAX_SENSOR_COEFF= 0.15\nADJ_JST_SENSOR_COEFF= ( 0.5, 0.02 )\nCONV_NUM_METHOD_FLOW= JST\nENTROPY_FIX_COEFF= 0.001\nTIME_DISCRE_FLOW= EULER_IMPLICIT\nRELAXATION_FACTOR_FLOW= 0.95\nCONV_NUM_METHOD_TURB= SCALAR_UPWIND\nTIME_DISCRE_TURB= EULER_IMPLICIT\nCFL_REDUCTION_TURB= 1.0\nRELAXATION_FACTOR_TURB= 0.95\nCONV_NUM_METHOD_ADJFLOW= JST\nTIME_DISCRE_ADJFLOW= EULER_IMPLICIT\nRELAXATION_FACTOR_ADJFLOW= 1.0\nCFL_REDUCTION_ADJFLOW= 0.8\nLIMIT_ADJFLOW= 1E6\nGEO_MARKER= ( AIRFOIL )\nGEO_DESCRIPTION= AIRFOIL\nGEO_MODE= FUNCTION\nCONV_CRITERIA= RESIDUAL\nRESIDUAL_REDUCTION= 9\nRESIDUAL_MINVAL= -12\nSTARTCONV_ITER= 10\nCAUCHY_ELEMS= 100\nCAUCHY_EPS= 1E-6\nCAUCHY_FUNC_FLOW= DRAG\nCAUCHY_FUNC_ADJFLOW= SENS_GEOMETRY\nMESH_FILENAME= mesh_RAE2822_turb.su2\nMESH_FORMAT= SU2\nMESH_OUT_FILENAME= mesh_out.su2\nSOLUTION_FLOW_FILENAME= solution_flow.dat\nSOLUTION_ADJ_FILENAME= solution_adj.dat\nOUTPUT_FORMAT= PARAVIEW\nCONV_FILENAME= history\nRESTART_FLOW_FILENAME= restart_flow.dat\nRESTART_ADJ_FILENAME= restart_adj.dat\nVOLUME_FLOW_FILENAME= flow\nVOLUME_ADJ_FILENAME= adjoin\nGRAD_OBJFUNC_FILENAME= of_grad.dat\nSURFACE_FLOW_FILENAME= surface_flow\nSURFACE_ADJ_FILENAME= surface_adjoint\nWRT_SOL_FREQ= 50.0\nWRT_CON_FREQ= 1\nDV_MARKER= ( AIRFOIL )\nOPT_OBJECTIVE= DRAG * 1.0\nOPT_CONSTRAINT= ( MOMENT_Z < 0.093 ) * 0.001; ( AIRFOIL_THICKNESS > 0.12 ) * 0.001\nOPT_GRADIENT_FACTOR= 1E-6\nOPT_RELAX_FACTOR= 1E2\nOPT_ITERATIONS= 100\nOPT_ACCURACY= 1E-10\nOPT_LINE_SEARCH_BOUND= 1E6\nOPT_BOUND_UPPER= 1E10\nOPT_BOUND_LOWER= -1E10\nOPT_COMBINE_OBJECTIVE= NO\nVALUE_OBJFUNC_FILENAME= of_eval.dat\nSIDESLIP_ANGLE= 0.0\nFREESTREAM_PRESSURE= 101325.0\nMULTIPOINT_WEIGHT= (1.0)\nMULTIPOINT_MACH_NUMBER= (0.729)\nMULTIPOINT_AOA= (2.31)\nMULTIPOINT_SIDESLIP_ANGLE= (0.0)\nMULTIPOINT_REYNOLDS_NUMBER= (6.5E6)\nMULTIPOINT_TARGET_CL= (0.724)\nMULTIPOINT_FREESTREAM_PRESSURE= (101325.0)\nMULTIPOINT_FREESTREAM_TEMPERATURE= (288.15)\nNUMBER_PART= 4\nNZONES= 1\nCONSOLE= CONCISE\nGRADIENT_METHOD= DISCRETE_ADJOINT'
                string_dvparam =''
                string_defdv =''
                string_dvold =''
	        for i in range(len(self.p)):
			string_dvparam+=' ( '+str(self.p[i,0])+', '+str(self.p[i,1])+' )'
			string_defdv+=' ( 1 , 1.0 | AIRFOIL  | '+str(self.p[i,0])+' , '+str(self.p[i,1])+' )'
			string_dvold+='0.0'
                        if (i<len(self.p)-1):
				string_dvparam+=' ;'
				string_defdv+=' ;'
				string_dvold+=', '
                
                #print string_dvparam
                #print string_defdv

		for i in range(self.n):
                	string_dvnew =''
                        #print self.B[i]
                        #print self.B[i].size
		        for k in range(self.B[i].size):
                		string_dvnew +=str(self.B[i,k])
                        	if (k<self.B[i].size-1):
					string_dvnew+=', '
			
			#print string_dvnew		
 	
			with open('config'+str(i)+'.cfg', "w" ) as c:
				c.write(config)
                	        c.write('\nDV_KIND= '+'HICKS_HENNE, '*(self.n-1)+'HICKS_HENNE')
				c.write('\nDV_PARAM= '+string_dvparam)
                        	c.write('\nDEFINITION_DV= '+string_defdv)
                        	c.write('\nDV_VALUE_OLD= '+string_dvold)
				c.write('\nFIN_DIFF_STEP= '+'1E-6')                       
				c.write('\nDV_VALUE_NEW= '+string_dvnew)
				c.write('\nDV_VALUE= '+string_dvnew)
                # part for saving the results inside a sub-directory.                
                description='The file(s) config_i.cfg contain a runnable configuration file for simple 2D simulation.\n The Boundary Conditions are specified inside the lines 39-40 of the method WriteDraft.'
                involved_outputs=list()
                for i in range(self.n):
                    file_ith = str('config'+str(i)+'.cfg') 
                    involved_outputs.append(file_ith)
                invoked_method = 'WriteDraft'
                collection = CollectResults(class_name=self.class_name, involved_outputs=involved_outputs, invoked_method=invoked_method, description=description)
                collection.Collect()

	def RunSU2(self):
                """ this method will run the simulations: DEF module will be called before the CFD one.
		"""		
		for i in range(self.n):
                        if (i == 0):
                            os.system('rm -r CONFIG*')
			os.mkdir('CONFIG'+str(i))
			os.system('mv config'+str(i)+'.cfg CONFIG'+str(i)+'/')
			os.system('cp mesh_RAE2822_turb.su2 CONFIG'+str(i)+'/')
			os.system('cd CONFIG'+str(i)+'; SU2_DEF config'+str(i)+'.cfg >& output_def.txt; mv mesh_out.su2 mesh_RAE2822_turb.su2; SU2_CFD config'+str(i)+'.cfg >& output_cfd.txt')
		
                # part for saving the results inside a sub-directory.                
                description='The sub-folder(s) CONFIG contain all the results of the serial reun of a DEF module and CFD.\n The Boundary Conditions are specified inside the lines 39-40 of the method WriteDraft, and are the same listed into the configuration file inside every CONFIG sub-directory.'
                involved_outputs=list()
                for i in range(self.n):
                    file_ith = str('-r CONFIG'+str(i)) 
                    involved_outputs.append(file_ith)
                    print 'the ith name is:', involved_outputs[i]
                invoked_method = 'RunSU2'
                collection = CollectResults(class_name=self.class_name, involved_outputs=involved_outputs, invoked_method=invoked_method, description=description)
                collection.Collect()



        def RunADJ(self):
                """ This method will run the simulatons with the adjoint operator
                """
                files = self.ScanFolder()
                if ('.dat' in files[0]):
                    os.system('rm *.dat')
                if ('.csv' in files[0]):
                    os.system('rm *.csv')
                if ('.vtk' in files[0]):
                    os.system('rm *.vtk')

                for i in range(self.n):
                    print '--------------------------------------------------'
                    print 'we are into the cycle '+str(i)
                    if (i == 0):
                        os.system('rm -r CONFIG*')
                    os.mkdir('CONFIG'+str(i))
                    # ire : instead of doing only a mv config, also a copy is done.
                    os.system('cp config'+str(i)+'.cfg CONFIG'+str(i)+'/automatica.cfg')
                    os.system('mv config'+str(i)+'.cfg CONFIG'+str(i))
                    os.system('cp mesh_RAE2822_turb.su2 CONFIG'+str(i)+'/')
                    # tiziano's 1st version:
                    #os.system('cp AUTO.sh CONFIG'+str(i)+'/')
                    #os.system('cd CONFIG'+str(i)+'; SU2_DEF config'+str(i)+'.cfg >& output_def.txt; mv mesh_out.su2 mesh_RAE2822_turb.su2; ./AUTO.sh > auto_out.txt')
                    # ire 2nd draft:
                    #os.system('cd CONFIG'+str(i)+'/') 
                    res_from_cfd = 'forces_breakdown.dat history.csv restart_flow.dat'
                    absolute = self.ScanFolder()
                    actual   = str(absolute[1])+'/CONFIG'+str(i)
                    if (i > 0):
                        os.system('SU2_DEF '+actual+'/config'+str(i)+'.cfg >& '+actual+'/output_def.txt')
                        os.system('mv mesh_out.su2 '+actual)
                        os.system('cp '+actual+'/mesh_out.su2 '+actual+'/mesh_RAE2822_turb.su2')
                    with open(actual+"/AUTO.sh", "w") as auto:
                            auto.write('#!/bin/bash')
                            auto.write('\n#------SCRIPT FOR AUTOMATIC DIRECT & ADJOINT CALCULATIONS AND GRADIENTS EVALUATION------#')
                            auto.write('\n#------CALCULATING DIRECT SOLUTION------#')
                            auto.write('\nSU2_CFD '+actual+'/automatica.cfg')
                            auto.write('\nmv '+res_from_cfd+' '+actual)
                            auto.write('\nSU2_SOL '+actual+'/automatica.cfg')
                            auto.write('\nmv flow.vtk '+actual)
                            auto.write('\n#------GENERATING NEW CONFIGURATION FILE FOR ADJOINT CALCULATIONS------#')
                            auto.write('\nsed -e "3,75s/DIRECT/DISCRETE\_ADJOINT/" '+actual+'/automatica.cfg > '+actual+'/ADJOINT.cfg')
                            auto.write('\n#------RENAMING DIRECT SOLUTION FILE TO ALLOW FOR RECOGNITION BY ADJOINT SOLVER------#')
                            auto.write('\ncp '+actual+'/restart_flow.dat '+actual+'/solution_flow.dat')
                            auto.write('\n#------CALCULATING ADJOINT SOLUTION------#')
                            auto.write('\nSU2_CFD_AD '+actual+'/ADJOINT.cfg > '+actual+'/output_adj.txt')
                            auto.write('\nmv surface_flow.csv surface_flow.vtk '+actual)
                            auto.write('\n#------GENERATING NEW CONFIGURATION FILE FOR SENSITIVITY PROJECTION------#')
                            auto.write('\nsed -e "74s/history\_adjoint/history\_dot\_ad/" '+actual+'/ADJOINT.cfg > '+actual+'/DOT_AD.cfg')
                            auto.write('\n#------RENAMING ADJOINT SOLUTION FILE TO ALLOW FOR RECOGNITION BY ALGORITHMIC DIFFERENCIATOR------#')
                            # 11/12/2018: the following lines give problems because the file restart_adj_cd.dat is not found.
                            #auto.write('\nmv restart_adj_cd.dat '+actual)
                            #auto.write('\ncp '+actual+'/restart_adj_cd.dat '+actual+'/solution_adj_cd.dat')
                            auto.write('\n#------GRADIENTS EVALUATION BY SURFACE SENSITIVITY PROJECTION INTO THE DESIGN SPACE------#')
                            auto.write('\nSU2_DOT_AD '+actual+'/DOT_AD.cfg')
                           
                    os.system('chmod +x '+actual+'/AUTO.sh')
                    os.system('chmod 755 '+actual+'/AUTO.sh')
                    os.system(actual+'/AUTO.sh > '+actual+'/auto_out.txt')

        def ReadRes(self):
		""" this method will read the results
		"""		
                cd=list()
		for i in range(self.n):
			results_file='CONFIG'+str(i)+'/forces_breakdown.dat'
			with open(results_file, 'r') as f:
				lines = f.read().split("\n")
                                check = False
				for k, line in enumerate(lines):
					if 'Surface name: AIRFOIL' in line:	
						check = True
						print 'the number of row for Surface Airfoil is:', k+1
					if ('Total CD' in line) & (check):
						print 'reading the Drag in line', k+1
						row_cd = k+1
						l = lines[k]
                                               	entries=l.split()
                                                cd.append(entries[4])
                for k in range(len(cd)):
                    cd[k] = float(cd[k])
                with open("cd_array.txt", "w") as grad:
                    grad.write('%s' % cd)
                gradients = np.zeros(len(cd)-1)
                for k in range(len(cd)-1):
                    cdi = (cd[k+1])
                    gradient = (cdi-cd[0])/self.B[k+1,k]
                    #print 'the gradient is:', gradient
                    gradients[k]= gradient
                with open("gradients_results.txt","w") as res:
                    res.write('%s' % gradients)
                        


		return cd, gradients
