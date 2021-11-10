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
            'folder'        : '/disk/zbox/ttamfa/dustevo/3DY1-42/',
            'nghost'        : 2,
            'xmax'          : 4*256,
            'ymax'          : 4*256,
            'ipart'         : 1,                 # WITH PARTICLES   = 1
            'igrowth'       : 0,                 # WITH DUST GROWTH = 1
            'nx_nodes'      : 1,
            'ny_nodes'      : 1,
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
            'sig_ref'       : 53.7587,
            'v_ref'         : 1010, 
            'value_label'   : 'DY1-42'}

arguments2 = arguments1.copy()
arguments2['folder']        = '/disk/zbox/ttamfa/dustevo/RY1-12/'
#arguments2['xmax']          = 4*256
#arguments2['ymax']          = 4*256
#arguments2['nx_nodes']      = 1
#arguments2['ny_nodes']      = 1 
arguments2['value_label']   = 'RY1-12'


arguments3 = arguments1.copy()
arguments3['folder']        = '/disk/zbox/ttamfa/dustevo/RY1-22/'
#arguments3['xmax']          = 1*256
#arguments3['ymax']          = 1*256
#arguments3['nx_nodes']      = 4
#arguments3['ny_nodes']      = 4 
arguments3['value_label']   = 'RY1-22'

arguments4 = arguments1.copy()
arguments4['folder']        = '/disk/zbox/ttamfa/dustevo/DY32/'
#arguments4['xmax']          = 1*256
#arguments4['ymax']          = 1*256
#arguments4['nx_nodes']      = 4
#arguments4['ny_nodes']      = 4 
arguments4['value_label']   = 'DY32'

arguments5 = arguments1.copy()
#arguments5['folder']        = '/Users/ttamfa/Desktop/TESTUNGROSSBI/TEST2/'
#arguments5['xmax']          = 4*256
#arguments5['ymax']          = 4*256
#arguments5['nx_nodes']      = 1
#arguments5['ny_nodes']      = 1 
#arguments5['value_label']   = 'blaefre5'


def loopover(loop, arguments1, arguments2, arguments3,arguments4,arguments5):
    data_stat1   = rf.read_stationary_data(arguments1)
    data1        = rf.read_single_timestep(loop,arguments1)
    
    data_stat2   = rf.read_stationary_data(arguments2)
    data2        = rf.read_single_timestep(loop,arguments2)
    
    data_stat3   = rf.read_stationary_data(arguments3)
    data3        = rf.read_single_timestep(loop,arguments3)

    #data_stat4   = rf.read_stationary_data(arguments4)
    #data4        = rf.read_single_timestep(loop,arguments4)

    #data_stat5   = rf.read_stationary_data(arguments5)
    #data5        = rf.read_single_timestep(loop,arguments5)


    rpf.plot_singlehist(loop,  
                    arguments1, 
                    data_stat1['r'], 
                    data1['rosp']*arguments1['sig_ref'], 
                    data1['rosg']*arguments1['sig_ref'], 
                    data1_dust_norm     =data_stat1['rospk']*arguments1['sig_ref'], 
                    data1_gas_norm      =data_stat1['rosk']*arguments1['sig_ref'],
                    plotting_arguments2 = arguments2, 
                    r2                  = data_stat2['r'], 
                    data2_dust          = data2['rosp']*arguments2['sig_ref'], 
                    data2_gas           = data2['rosg']*arguments2['sig_ref'],
                    data2_dust_norm     = data_stat2['rospk']*arguments2['sig_ref'], 
                    data2_gas_norm      = data_stat2['rosk']*arguments2['sig_ref'], 
                    plotting_arguments3 = arguments3, 
                    r3                  = data_stat3['r'], 
                    data3_dust          = data3['rosp']*arguments3['sig_ref'], 
                    data3_gas           = data3['rosg']*arguments3['sig_ref'], 
                    data3_dust_norm     = data_stat3['rospk']*arguments3['sig_ref'], 
                    data3_gas_norm      = data_stat3['rosk']*arguments3['sig_ref'], 
                    plotting_arguments4 = arguments4, 
                    r4                  = [],#data_stat4['r'], 
                    data4_dust          = [],#data4['rosp']*arguments4['sig_ref'], 
                    data4_gas           = [],#data4['partsize'],#*arguments4['sig_ref'], 
                    data4_dust_norm=[], 
                    data4_gas_norm=[], 
                    plotting_arguments5 = arguments5, 
                    r5=[], 
                    data5_dust=[], 
                    data5_gas=[], 
                    data5_dust_norm=[], 
                    data5_gas_norm=[], 
                    planetpos=[], plotrangedust=[], plotrangegas=[])


loopover(401,arguments1, arguments2, arguments3, arguments4, arguments5)

#num_cores = multiprocessing.cpu_count()
#Parallel(n_jobs=num_cores)(delayed(loopover)(i, arguments) for i in range(1,4))


