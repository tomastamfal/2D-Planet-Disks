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


#from joblib import Parallel, delayed
#import multiprocessing



# =====================================================
#  Definitions
# =====================================================

arguments1 = {
            'folder'        : '/Users/ttamfa/Desktop/DY32/',
            'nghost'        : 2,
            'xmax'          : 1*256,
            'ymax'          : 1*256,
            'ipart'         : 1,                 # WITH PARTICLES   = 1
            'igrowth'       : 1,                 # WITH DUST GROWTH = 1
            'nx_nodes'      : 4,
            'ny_nodes'      : 4,
            'R_ref'         : 10,                 # IN AU
            'planet'        : 'False',
            'bg_color'      : 'white',
            'fg_color'      : 'black',
            'label_color'   : 'black',
            'colormap1'     : 'hot',
            'colormap2'     : 'plasma',
            'colormap3'     : 'viridis',
            'name'          : 'compare',
            'fontsize'	    : 12,	
	    'line_width'    : 5.0,
            'sig_ref'       : 53.7587,
            'v_ref'         : 1010}

arguments2 = {
            'folder'        : '/Users/ttamfa/Desktop/RY32/',
            'nghost'        : 2,
            'xmax'          : 1*256,
            'ymax'          : 1*256,
            'ipart'         : 1,                 # WITH PARTICLES   = 1
            'igrowth'       : 0,                 # WITH DUST GROWTH = 1
            'nx_nodes'      : 4,
            'ny_nodes'      : 4,
            'R_ref'         : 10,                 # IN AU
            'a'             : 3,
            'planet'        : 'False',
            'bg_color'      : 'white',
            'fg_color'      : 'black',
            'label_color'   : 'black',
            'colormap1'     : 'hot',
            'colormap2'     : 'plasma',
            'colormap3'     : 'viridis',
            'name'          : 'compare',
            'fontsize'	    : 12,	
	    'line_width'    : 5.0,
            'sig_ref'       : 53.7587,
            'v_ref'         : 1010}



def loopover(loop,arguments1, arguments2):
    data_stat1   = rf.read_stationary_data(arguments1)
    data1        = rf.read_single_timestep(loop,arguments1)
    
    data_stat2   = rf.read_stationary_data(arguments2)
    data2        = rf.read_single_timestep(loop,arguments2)
 

    stokes1 = np.pi*0.5*data1['partsize']/(data1['rosg']*arguments1['sig_ref'])
    stokes2 = np.pi*0.5*arguments2['a']/(data2['rosg']*arguments2['sig_ref'])

    rpf.plot_compare_2D(loop, arguments1, 
		    	data_stat1['r'], 
                        data_stat1['theta'], 
                        data1['rosp']/data_stat1['rospk'], 
                        data1['rosg']/data_stat1['rosk'],
                        stokes1,
                        data1['partsize'],
                        data1['regime'],
                        arguments2, 
                        data_stat2['r'], 
                        data_stat2['theta'], 
                        data2['rosp']/data_stat2['rospk'], 
                        data2['rosg']/data_stat2['rosk'],
                        stokes2,
                        plotrange=[1e-8, 1e4, 1e-1, 1e1, 1e-6, 1e-1])

#plotrange=[0.5, 3, 1e-1, 1e1, 1e-6, 1e-1]
for kk in range(201,301):
    loopover(kk,arguments1,arguments2)

#num_cores = multiprocessing.cpu_count()
#Parallel(n_jobs=num_cores)(delayed(loopover)(i, arguments) for i in range(1,4))


