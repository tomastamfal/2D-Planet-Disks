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
            'folder'        : '/disk/zbox/ttamfa/dustevo/ALL_ZHU/ZHU',
            'nghost'        : 2,
            'xmax'          : 1*282,
            'ymax'          : 1*1024,
            'ipart'         : 1,                 # WITH PARTICLES   = 1
            'igrowth'       : 0,                 # WITH DUST GROWTH = 1
            'nx_nodes'      : 1,
            'ny_nodes'      : 1,
            'R_ref'         : 1.,                 # IN AU
            'rin'           : 0.5,
            'rout'          : 3,
            'planet'        : 'True',
            'bg_color'      : 'white',
            'fg_color'      : 'black',
            'label_color'   : 'black',
            'fontsize'	    : 14,	
	    'line_width'    : 1.0,
            'name'          : 'comparelineplot',
            'value_name'    : r'log($\Sigma(r,t)$ / $\Sigma_0(r=1, t=0)$)',
            'sig_ref'       : 178,
            'v_ref'         : 1010, 
            'value_label'   : '282'}


arguments2 = arguments1.copy()
arguments2['folder']        = '/disk/zbox/ttamfa/dustevo/DZHU'
arguments2['xmax']          = 1*282
#arguments2['ymax']          = 4*256
#arguments2['nx_nodes']      = 1
#arguments2['ny_nodes']      = 1 
arguments2['value_label']   = '282 (D)'


arguments3 = arguments1.copy()
arguments3['folder']        = '/disk/zbox/ttamfa/dustevo/ALL_ZHU/ZHU'
#arguments3['xmax']          = 1*2048
#arguments3['ymax']          = 1*256
#arguments3['nx_nodes']      = 4
#arguments3['ny_nodes']      = 4 
arguments3['value_label']   = 'Initial'

arguments4 = arguments1.copy()
#arguments4['folder']        = '/disk/zbox/ttamfa/dustevo/ZHU_DUSTEVO/'
#arguments4['xmax']          = 1*2084
#arguments4['ymax']          = 1*256
#arguments4['nx_nodes']      = 4
#arguments4['ny_nodes']      = 4 
arguments4['value_label']   = 'Zhu et al.'

arguments5 = arguments1.copy()
#arguments5['folder']        = '/disk/zbox/ttamfa/dustevo/ZHU_DUSTEVO2'
#arguments5['xmax']          = 1*282
#arguments5['ymax']          = 4*256
#arguments5['nx_nodes']      = 1
#arguments5['ny_nodes']      = 1 
arguments5['value_label']   = '282 (D) WRONG'


def loopover(loop, arguments1, arguments2=[], arguments3=[],arguments4=[],arguments5=[]):
    R_ref        = arguments1['R_ref']         
    
    data_stat1   = rf.read_stationary_data(arguments1)
    data1        = rf.read_single_timestep(loop,arguments1)
    datanorm     = rf.read_single_timestep(1,arguments1)
    #datanorm     = data_stat1

    data_index  = rf.find_index(data_stat1['r'][:,1], 1.*R_ref)
    print data_index
    print data_stat1['r'][data_index,:]
    print data_stat1['r'][:,1]
    data_0      = datanorm['rosp']
    data_tmp    =    data1['rosp']
    data1_norm  = 10000000000000000000000000000000000.0*data_tmp
    for j in range(np.shape(data_tmp)[1]):
        data1_norm[:,j] = data_tmp[:,j]/data_0[data_index,j]


    data_stat2   = rf.read_stationary_data(arguments2)
    data2        = rf.read_single_timestep(loop,arguments2)
    datanorm     = rf.read_single_timestep(1,arguments2)
    #datanorm     = data_stat2

    data_index  = rf.find_index(data_stat2['r'][:,1], 1.*R_ref)
    data_0      = datanorm['rosp']
    data_tmp    =    data2['rosp']
    data2_norm  = 0.0*data_tmp
    for j in range(np.shape(data_tmp)[1]):
        data2_norm[:,j] = data_tmp[:,j]/data_0[data_index,j]

    
    data_stat3   = rf.read_stationary_data(arguments3)
    data3        = rf.read_single_timestep(1,arguments3) 
    datanorm     = data3 

    data_index  = rf.find_index(data_stat3['r'][:,1], 1.*R_ref)
    data_0      = datanorm['rosp']
    data_tmp    =    data3['rosp']
    data3_norm  = 0.0*data_tmp
    for j in range(np.shape(data_tmp)[1]):
        data3_norm[:,j] = data_tmp[:,j]/data_0[data_index,j]

    tmp     = np.loadtxt('plot2_part_b.txt',  delimiter=',')
    rzhu    = tmp[:,0]*R_ref
    datazhu = tmp[:,1]
    #data_stat4   = rf.read_stationary_data(arguments4)
    #data4        = rf.read_single_timestep(loop,arguments4)
    #datanorm     = rf.read_single_timestep(1,arguments4)

    #data_index  = rf.find_index(data_stat4['r'][:,1], 1)
    #data_0      = datanorm['rosp']
    #data_tmp    =    data4['rosp']
    #data4_norm  = 0.0*data_tmp
    #for j in range(np.shape(data_tmp)[1]):
    #    data4_norm[:,j] = data_tmp[:,j]/data_0[data_index,j]


    #data_stat5   = rf.read_stationary_data(arguments5)
    #data5        = rf.read_single_timestep(loop,arguments5)
    #datanorm     = rf.read_single_timestep(1,arguments5)

    #data_index  = rf.find_index(data_stat5['r'][:,1], 1)
    #data_0      = datanorm['rosp']
    #data_tmp    =    data5['rosp']
    #data5_norm  = 0.0*data_tmp
    #for j in range(np.shape(data_tmp)[1]):
    #    data5_norm[:,j] = data_tmp[:,j]/data_0[data_index,j]




    rpf.plot_singlelineplotzhu(loop,  
                    arguments1, 
                    data_stat1['r'], 
                    np.log10(data1_norm), 
                    plotting_arguments2 = arguments2, 
                    r2                  = data_stat2['r'], 
                    data2               = np.log10(data2_norm), 
                    plotting_arguments3 = arguments3, 
                    r3                  = data_stat3['r'], 
                    data3               = np.log10(data3_norm), 
                    plotting_arguments4 = arguments4, 
                    r4                  = rzhu, 
                    data4               = datazhu, 
                    plotting_arguments5 = arguments5, 
                    r5=[],#data_stat5['r'], 
                    data5=[],#np.log10(data5_norm), 
                    planetpos=[], plotrange=[0.6, 1.5, -0.6, 1.2])


loopover(201,arguments1, arguments2, arguments3, arguments4, arguments5)

#num_cores = multiprocessing.cpu_count()
#Parallel(n_jobs=num_cores)(delayed(loopover)(i, arguments) for i in range(1,4))


