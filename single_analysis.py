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
import numpy as np
import rossbi_functions as rf
import rossbi_plottingfunctions as rpf


#from joblib import Parallel, delayed
#import multiprocessing



# =====================================================
#  Definitions
# =====================================================

arguments1 = {
            'folder'        : '/Users/ttamfa/Desktop/GEO_NO_YES/',
            'nghost'        : 2,
            'xmax'          : 1*256,
            'ymax'          : 1*256,
            'ipart'         : 1,                 # WITH PARTICLES   = 1
            'igrowth'       : 0,                 # WITH DUST GROWTH = 1
            'nx_nodes'      : 8,
            'ny_nodes'      : 4,
            'R_ref'         : 1,                 # IN AU
            'rin'           : 0.3,
            'rout'          : 10,
            'bg_color'      : 'white',
            'fg_color'      : 'black',
            'label_color'   : 'black',
            'fontsize'	    : 10,	
	    'line_width'    : 1,
            'name'          : 'bdsjhbfhsdbfhbsdi',
            'sig_ref'       : 178,
            'v_ref'         : 1010, 
            'value_label'   : 'TEST2'}


def loopover(loop, arguments1):
    data_stat1   = rf.read_stationary_data(arguments1)
    data1        = rf.read_single_timestep(loop,arguments1)
    datanorm     = rf.read_single_timestep(588:q
            ,arguments1)


    rpf.single_analysis(loop,  
                    data_stat1['r'], 
                    data1['rosp']*arguments1['sig_ref'], 
                    data1['rosg']*arguments1['sig_ref'], arguments1,
                    data1_dust_norm = datanorm['rosp']*arguments1['sig_ref'], 
                    data1_gas_norm  = datanorm['rosg']*arguments1['sig_ref'],)


loopover(588,arguments1)

#num_cores = multiprocessing.cpu_count()
#Parallel(n_jobs=num_cores)(delayed(loopover)(i, arguments) for i in range(1,4))


