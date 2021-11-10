# =====================================================
#   The variable names are  
# =====================================================
#   ug_glob,    up_glob 		are radial gas or particle velocity
#   vg_glob,    vp_glob 		are azimuthal gas or particle velocity
#   rosg_glob,  rosp_glob 	        are gas or particle density
#   pg_glob,    pg_glob 		is gas pressure
#   vort_glob,  vortg_glob 	        is gas vorticity
#
#   vk_glob,    vpk_glob 		are gas or particle background azim. vel.
#   rosk_glob,  rospk_glob 	        are gas or particle backgr. density
#   pk_glob,    pk_glob 		is gas background pressure
#   vortk_glob, vortk_glob 	        is gas background vorticity
#
#   partsize 
#   stokes
#   regime
#
# =====================================================
import rossbi_functions as rf
import rossbi_plottingfunctions as rpf
import numpy as np

#from joblib import Parallel, delayed
#import multiprocessing

#  Definitions
# =====================================================

arguments1 = {
            'folder'        : '/Users/ttamfa/Desktop/DY1-42/',
            'nghost'        : 2,
            'xmax'          : 1*256,
            'ymax'          : 1*256,
            'ipart'         : 1,                 # WITH PARTICLES   = 1
            'igrowth'       : 0,                 # WITH DUST GROWTH = 1
            'nx_nodes'      : 4,
            'ny_nodes'      : 4,
            'R_ref'         : 10,                 # IN AU
            'rin'           : 0.5,
            'rout'          : 3,
            'planet'        : 'False',
            'bg_color'      : 'white',
            'fg_color'      : 'black',
            'label_color'   : 'black',
            'fontsize'	    : 10,	
	    'line_width'    : 1.0,
            'name'          : 'new',
            'dust_gas'      : 1e-2,
            'a0'            : 1e-4, 
            'sig_ref'       : 53.7587,
            'v_ref'         : 10**(-0.5) * 2.97396e6, 
            'value_label'   : 'DY32',
            'rho_sol'       : 1.0}


def loopover(loop, arguments1):
    data_stat1   = rf.read_stationary_data(arguments1)
    data1        = rf.read_single_timestep(loop,arguments1)
    
    time = np.loadtxt('time.data')
    time = time[loop]
    gamma = 1.4
    rpf.partsize_analysis(loop,  
                    arguments1, 
                    data_stat1['r'],
                    data_stat1['theta'],
                    data1['rosp'], 
                    data1['rosg'], 
                    data1['pg'],
                    np.sqrt(gamma*data1['pg']/(data1['rosg'])),
                    time,
                    planetpos=[],
                    plotrange=[])
                    


loopover(401,arguments1)

#num_cores = multiprocessing.cpu_count()
#Parallel(n_jobs=num_cores)(delayed(loopover)(i, arguments) for i in range(1,4))


