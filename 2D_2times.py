# =====================================================
#   The variable names are  
# =====================================================
#   ug_glob,    up_glob 		are radial gas or particle velocity
#   vg_glob,    vp_glob 		are azimuthal gas or particle velocity
#   rosg_glob,  rosp_glob 	are gas or particle density
#   pg_glob,    pg_glob 		is gas pressure
#   vort_glob,  vortg_glob 	is gas vorticity
#
#   vk_glob,    vpk_glob 		are gas or particle background azim. vel.
#   rosk_glob,  rospk_glob 	are gas or particle backgr. density
#   pk_glob,    pk_glob 		is gas background pressure
#   vortk_glob, vortk_glob 	is gas background vorticity
#
#   partsize 
#   stokes
#   regime
#
# =====================================================
import numpy as np
import rossbi_functions as rf
import rossbi_plottingfunctions as rpf


from joblib import Parallel, delayed
import multiprocessing



# =====================================================
#  Definitions
# =====================================================

arguments = {
            'folder'        : '/disk/zbox/ttamfa/dustevo/RY32/',
            'nghost'        : 2,
            'xmax'          : 4*256,
            'ymax'          : 4*256,
            'ipart'         : 1,                 # WITH PARTICLES   = 1
            'igrowth'       : 0,                 # WITH DUST GROWTH = 1
            'nx_nodes'      : 1,
            'ny_nodes'      : 1,
            'R_ref'         : 10,  
            'planet'        : 1,
            'bg_color'      : 'white',
            'fg_color'      : 'black',
            'label_color'   : 'black',
            'colormap1'     : 'hot',
            'colormap2'     : 'viridis_r',
            'colormap3'     : 'viridis_r',
            'fontsize'      : 8, 
	    'line_width'    : 0.5,
            'name'          : '2D_2times',
            'sig_ref'       : 53.7587,
            'v_ref'         : 1010}



def loopover(loop,arguments):
    data_stat   = rf.read_stationary_data(arguments)
    data        = rf.read_single_timestep(loop,arguments)
    
    rpf.plot_single2D_2times(loop, arguments, 
		    		data_stat['r'], 
                                data_stat['theta'], 
                                data['rosp']/data_stat['rospk'], 
                                data['rosg']/data_stat['rosk'])

loopover(400,arguments)

#num_cores = multiprocessing.cpu_count()
#Parallel(n_jobs=num_cores)(delayed(loopover)(i, arguments) for i in range(1,4))


