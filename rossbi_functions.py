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

import os 
import numpy as np

# =====================================================
#  Definitions
# =====================================================

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def read_stationary_data(arguments):

    nghost          = arguments['nghost']
    xmax            = arguments['xmax'] + 2*nghost
    ymax            = arguments['ymax'] + 2*nghost
    ipart           = arguments['ipart']                 
    igrowth         = arguments['igrowth']                 
    nx_nodes        = arguments['nx_nodes']
    ny_nodes        = arguments['ny_nodes']
    folder          = arguments['folder']

    
    os.chdir(folder)
    if (arguments['planet'] == 'True'):
        planet  = np.loadtxt('planet.data')
    else:
        planet  = np.zeros(10)
    
    # Create global variable: for the total grid
    xmax_glob   = nx_nodes*(xmax-2*nghost) + 2*nghost
    ymax_glob   = ny_nodes*(ymax-2*nghost) + 2*nghost


    nx = xmax-2*nghost
    ny = ymax-2*nghost

    r_glob      = np.zeros((xmax_glob, ymax_glob))
    theta_glob  = np.zeros((ymax_glob, 1))
    dr_glob     = np.zeros((xmax_glob, 1))
    dtheta_glob = 0


    #BACKGROUND GAS
    rosk_glob   = np.zeros((xmax_glob, ymax_glob))
    pk_glob     = np.zeros((xmax_glob, ymax_glob))
    uk_glob     = np.zeros((xmax_glob, ymax_glob))
    vk_glob     = np.zeros((xmax_glob, ymax_glob))
    vortk_glob  = np.zeros((xmax_glob, ymax_glob))

    # PARTICLES
    if ipart != 0:
        #BACKGROUND
        rospk_glob = np.zeros((xmax_glob, ymax_glob))
        upk_glob   = np.zeros((xmax_glob, ymax_glob))
        vpk_glob   = np.zeros((xmax_glob, ymax_glob))
   	vortpk_glob= np.zeros((xmax_glob, ymax_glob))


    #START WITH THE MPI STUFF:
    for i_node in range(1,nx_nodes+1):
        i_node_str = '%3.3d'% i_node

        for j_node in range(1,ny_nodes+1):
            icount=1
            j_node_str = '%3.3d'  %j_node

            grid_pos   = '.%s_%s' %(i_node_str, j_node_str)

            

            time=np.loadtxt('time.data')

            r        = np.zeros((xmax,ymax))
            theta    = np.zeros((ymax,1))
            dr       = np.zeros((xmax,1))
            dtheta   = 0
            
            rosk     = np.zeros((xmax,ymax))
            uk       = np.zeros((xmax,ymax))
            vk       = np.zeros((xmax,ymax))
            vortk    = np.zeros((xmax,ymax))
            pk       = np.zeros((xmax,ymax))
            
            
            # Import of r
            with open('r.data%s' %grid_pos, 'rb') as f:
                f.read(4)
                Y = np.fromfile(f, dtype=np.float64)
                Y = np.append(Y, [0])
                r[:,0] = Y[0::2]    
	    r = (np.ones((ymax,1))*r[:,0]).transpose()

            # Import of theta
            with open('theta.data%s' %grid_pos, 'rb') as f:
                f.read(4)
                Y = np.fromfile(f, dtype=np.float64)
                Y = np.append(Y, [0])
                theta[:,0] = Y[0::2]

            dtheta = theta[1]-theta[0]

            # Stationnary variables
            with open('Stationnary_gas.data%s' %grid_pos, 'rb') as f:
                f.read(4)
                Y = np.fromfile(f, dtype=np.float64)
                Y = np.append(Y, [0])

                uk[:,0]    = Y[0:5*(xmax):5]
                vk[:,0]    = Y[1:5*(xmax):5]
                pk[:,0]    = Y[2:5*(xmax):5]
                rosk[:,0]  = Y[3:5*(xmax):5]

            uk   = (np.ones((ymax,1))*uk[:,0]).transpose()
            vk   = (np.ones((ymax,1))*vk[:,0]).transpose()
            pk   = (np.ones((ymax,1))*pk[:,0]).transpose()
            rosk = (np.ones((ymax,1))*rosk[:,0]).transpose()


            beta_omega=-3/2
            for i in range(xmax-1):
                dr[i]=r[i+1,1]-r[i,1]
            dr[xmax-1]=dr[xmax-2]


            for i in range(1,xmax-1):
                vortk[i,:] = 1/r[i,:]*((r[i+1,:]*vk[i+1,:]-r[i-1,:]*vk[i-1,:])/(dr[i]+dr[i-1]))
               

            r_glob[np.ix_((i_node-1)*nx + np.array(range(xmax)), (j_node-1)*ny + np.array(range(ymax)))] = r
	    theta_glob[np.ix_((j_node-1)*ny + np.array(range(ymax)))] = theta
	    dr_glob[np.ix_((i_node-1)*nx + np.array(range(xmax)))] = dr
	    vk_glob[np.ix_((i_node-1)*nx + np.array(range(xmax)), (j_node-1)*ny + np.array(range(ymax)))]            = vk
	    rosk_glob[np.ix_((i_node-1)*nx + np.array(range(xmax)), (j_node-1)*ny + np.array(range(ymax)))]          = rosk
            pk_glob[np.ix_((i_node-1)*nx + np.array(range(xmax)), (j_node-1)*ny + np.array(range(ymax)))]            = pk
	    vortk_glob[np.ix_((i_node-1)*nx + np.array(range(1,xmax-1)), (j_node-1)*ny + np.array(range(1,ymax-1)))] = vortk[np.ix_(np.array(range(1,xmax-1)), np.array(range(1,ymax-1)))]
         
	    if ipart !=0:
                rospk    = np.zeros((xmax,ymax))
                upk      = np.zeros((xmax,ymax))
                vpk      = np.zeros((xmax,ymax))

                with open('Stationnary_part.data%s' %grid_pos, 'rb') as f:
                    f.read(4)
                    Y = np.fromfile(f, dtype=np.float64)
                    Y = np.append(Y, [0])
                    upk[:,0]    = Y[0:4*(xmax):4]
                    vpk[:,0]    = Y[1:4*(xmax):4]
                    rospk[:,0]  = Y[2:4*(xmax):4]


                vortpk  = vortk
                upk     = (np.ones((ymax,1))*upk[:,0]).transpose()
                vpk     = (np.ones((ymax,1))*vpk[:,0]).transpose()
                rospk   = (np.ones((ymax,1))*rospk[:,0]).transpose()
                

                vpk_glob[np.ix_((i_node-1)*nx + np.array(range(xmax)), (j_node-1)*ny + np.array(range(ymax)))]            = vpk
                rospk_glob[np.ix_((i_node-1)*nx + np.array(range(xmax)), (j_node-1)*ny + np.array(range(ymax)))]          = rospk
                vortpk_glob[np.ix_((i_node-1)*nx + np.array(range(1,xmax-1)), (j_node-1)*ny + np.array(range(1,ymax-1)))] = vortpk[np.ix_(np.array(range(1,xmax-1)), np.array(range(1,ymax-1)))]

    if (ipart!=0):
	return {'r'	    :r_glob,
		'dr'        :dr,		
		'theta'     :theta_glob,		
		'dtheta'    :dtheta,
		'vk'        :vk_glob,
                'rosk'      :rosk_glob,
                'pk'        :pk_glob,
                'vortk'     :vortk_glob,
                'vpk'       :vpk_glob,
                'rospk'     :rospk_glob,
                'vortpk'    :vortpk_glob,
                'planet'    :planet}


    if (ipart==0 ):
        return {'r'	    :r_glob,
		'dr'        :dr,		
		'theta'     :theta_glob,		
		'dtheta'    :dtheta,
		'vk'        :vk_glob,
                'rosk'      :rosk_glob,
                'pk'        :pk_glob,
                'vortk'     :vortk_glob, 
                'planet'    :planet}



def read_single_timestep(outerloop, arguments):
    nghost          = arguments['nghost']
    xmax            = arguments['xmax'] + 2*nghost
    ymax            = arguments['ymax'] + 2*nghost
    ipart           = arguments['ipart']                 
    igrowth         = arguments['igrowth']                 
    nx_nodes        = arguments['nx_nodes']
    ny_nodes        = arguments['ny_nodes']
    folder          = arguments['folder']
    
    os.chdir(folder)
    
    number          = '%04d' %outerloop    # Time index of the files

    # Create global variable: for the total grid
    xmax_glob   = nx_nodes*(xmax-2*nghost) + 2*nghost
    ymax_glob   = ny_nodes*(ymax-2*nghost) + 2*nghost


    #NORMAL GAS
    ug_glob     = np.zeros((xmax_glob, ymax_glob))
    vg_glob     = np.zeros((xmax_glob, ymax_glob))
    pg_glob     = np.zeros((xmax_glob, ymax_glob))
    rosg_glob   = np.zeros((xmax_glob, ymax_glob))
    vortg_glob  = np.zeros((xmax_glob, ymax_glob))
    divvg_glob  = np.zeros((xmax_glob, ymax_glob))

    omegak_glob = np.zeros((xmax_glob, ymax_glob))


    # PARTICLES
    if ipart != 0:
        #NORMAL
        up_glob     = np.zeros((xmax_glob, ymax_glob))
        vp_glob     = np.zeros((xmax_glob, ymax_glob))
        rosp_glob   = np.zeros((xmax_glob, ymax_glob))
        vortp_glob  = np.zeros((xmax_glob, ymax_glob)) 
        vortpk_glob = np.zeros((xmax_glob, ymax_glob)) 

        #GROWTH
        if igrowth != 0:
            regime_glob     = np.zeros((xmax_glob, ymax_glob))
            partsize_glob   = np.zeros((xmax_glob, ymax_glob))
            stokes_glob     = np.zeros((xmax_glob, ymax_glob))
            
          

    nx = xmax-2*nghost
    ny = ymax-2*nghost


    #START WITH THE MPI STUFF:
    for i_node in range(1,nx_nodes+1):
        i_node_str = '%3.3d'% i_node

        for j_node in range(1,ny_nodes+1):
            icount=1
            j_node_str = '%3.3d'  %j_node

            grid_pos   = '.%s_%s' %(i_node_str, j_node_str)
            file_gas   = '%s_gas.data%s' %(number,grid_pos)

            
            ug       = np.zeros((xmax,ymax))
            vg       = np.zeros((xmax,ymax))
            pg       = np.zeros((xmax,ymax))
            rosg     = np.zeros((xmax,ymax))
            vortg    = np.zeros((xmax,ymax))
            divvg    = np.zeros((xmax,ymax))
            


            with open('%s' %file_gas, 'rb') as f:
                f.read(4)
                Y = np.fromfile(f, dtype=np.float64)
                Y = np.append(Y, [0])

                ug    = np.reshape((Y[0:5*(ymax*xmax):5]), (ymax,xmax), order="F").transpose()
                vg    = np.reshape((Y[1:5*(ymax*xmax):5]), (ymax,xmax), order="F").transpose()
                pg    = np.reshape((Y[2:5*(ymax*xmax):5]), (ymax,xmax), order="F").transpose()
                rosg  = np.reshape((Y[3:5*(ymax*xmax):5]), (ymax,xmax), order="F").transpose()



            for i in range(1,xmax-1):
                for j in range(1,ymax-1):
                    #vortg[i,j] = 1/r[i,j]*(((r[i+1,j]*vg[i+1,j])-(r[i-1,j]*vg[i-1,j]))/(dr[i]+dr[i-1]) - ((ug[i,j+1]-ug[i,j-1])/(2*dtheta)))
	            vortg[i,j] = 1		

            for i in range(1,xmax-1):
                for j in range(1,ymax-1):
                    #divvg[i,j] = 1/r[i,j]*(((r[i+1,j]*ug[i+1,j])-(r[i-1,j]*ug[i-1,j]))/(dr[i]+dr[i-1]) - ((vg[i,j+1]-vg[i,j-1])/(2*dtheta)))
		    divvg[i,j] = 1	
            
	    
	    ug_glob[np.ix_((i_node-1)*nx + np.array(range(xmax)), (j_node-1)*ny + np.array(range(ymax)))]            = ug
            vg_glob[np.ix_((i_node-1)*nx + np.array(range(xmax)), (j_node-1)*ny + np.array(range(ymax)))]            = vg
            rosg_glob[np.ix_((i_node-1)*nx + np.array(range(xmax)), (j_node-1)*ny + np.array(range(ymax)))]          = rosg
            pg_glob[np.ix_((i_node-1)*nx + np.array(range(xmax)), (j_node-1)*ny + np.array(range(ymax)))]            = pg
            vortg_glob[np.ix_((i_node-1)*nx + np.array(range(1,xmax-1)), (j_node-1)*ny + np.array(range(1,ymax-1)))] = vortg[np.ix_(np.array(range(1,xmax-1)), np.array(range(1,ymax-1)))]
            divvg_glob[np.ix_((i_node-1)*nx + np.array(range(1,xmax-1)), (j_node-1)*ny + np.array(range(1,ymax-1)))] = divvg[np.ix_(np.array(range(1,xmax-1)), np.array(range(1,ymax-1)))]



            if ipart != 0:
                file_part='%s_part.data%s' %(number,grid_pos)
     
                if igrowth != 0:
                    file_growth='%s_growth.data%s' %(number,grid_pos)

                up       = np.zeros((xmax,ymax))
                vp       = np.zeros((xmax,ymax))
                rosp     = np.zeros((xmax,ymax))
                vortp 	 = np.zeros((xmax,ymax))
		divvp    = np.zeros((xmax,ymax))

                with open('%s' %file_part, 'rb') as f:
                    f.read(4)
                    Y = np.fromfile(f, dtype=np.float64)
                    Y = np.append(Y, [0])

                up   = np.reshape( (Y[0:4*(ymax*xmax):4]), (ymax,xmax), order="F").transpose()
                vp   = np.reshape( (Y[1:4*(ymax*xmax):4]), (ymax,xmax), order="F").transpose()
                rosp = np.reshape( (Y[2:4*(ymax*xmax):4]), (ymax,xmax), order="F").transpose()

                for i in range(1,xmax-1):
                    for j in range(1,ymax-1):
                        #vortp[i,j] = 1/r[i,j]*(((r[i+1,j]*vp[i+1,j]-r[i-1,j]*vp[i-1,j]))/(dr[i]+dr[i-1]) - ((up[i,j+1]-up[i,j-1])/(2*dtheta)))
			vortp[i,j] = 1

                for i in range(1,xmax-1):
                    for j in range(1,ymax-1):
                        #divvp[i,j] = 1/r[i,j]*(((r[i+1,j]*up[i+1,j]-r[i-1,j]*up[i-1,j]))/(dr[i]+dr[i-1]) - ((vp[i,j+1]-vp[i,j-1])/(2*dtheta)))
               		divvp[i,j] = 1 
                
                
                
                up_glob[np.ix_((i_node-1)*nx + np.array(range(xmax)), (j_node-1)*ny + np.array(range(ymax)))]             = up
                vp_glob[np.ix_((i_node-1)*nx + np.array(range(xmax)), (j_node-1)*ny + np.array(range(ymax)))]             = vp
                rosp_glob[np.ix_((i_node-1)*nx + np.array(range(xmax)), (j_node-1)*ny + np.array(range(ymax)))]           = rosp
                vortp_glob[np.ix_((i_node-1)*nx + np.array(range(1,xmax-1)), (j_node-1)*ny + np.array(range(1,ymax-1)))]  = vortp[np.ix_(np.array(range(1,xmax-1)), np.array(range(1,ymax-1)))]


                if igrowth != 0:
                    regime   = np.zeros((xmax,ymax))
                    partsize = np.zeros((xmax,ymax))
                    stokes   = np.zeros((xmax,ymax))

                    with open('%s' %file_growth, 'rb') as f:
                        f.read(4)
                        Y = np.fromfile(f, dtype=np.float64)
                        Y = np.append(Y, [0])

                    partsize   = np.reshape( (Y[0:4*(ymax*xmax):4]), (ymax,xmax), order="F").transpose()
                    regime     = np.reshape( (Y[1:4*(ymax*xmax):4]), (ymax,xmax), order="F").transpose()
                    stokes     = np.reshape( (Y[2:4*(ymax*xmax):4]), (ymax,xmax), order="F").transpose()
		    
                    regime_glob[np.ix_((i_node-1)*nx + np.array(range(nghost,nx+nghost)), (j_node-1)*ny + np.array(range(nghost,ny+nghost)))] 		= regime[nghost:-nghost, nghost:-nghost] 
		    partsize_glob[np.ix_((i_node-1)*nx + np.array(range(nghost,nx+nghost)), (j_node-1)*ny + np.array(range(nghost,ny+nghost)))] 	= partsize[nghost:-nghost, nghost:-nghost] 
		    stokes_glob[np.ix_((i_node-1)*nx + np.array(range(nghost,nx+nghost)), (j_node-1)*ny + np.array(range(nghost,ny+nghost)))]   	= stokes[nghost:-nghost, nghost:-nghost]
			    

    if (ipart!=0 and igrowth !=0):
        return {'ug'        :ug_glob, 
                'vg'        :vg_glob,
                'rosg'      :rosg_glob,
                'pg'        :pg_glob,
                'vortg'     :vortg_glob,
                'divvg'     :divvg_glob,
                'up'        :up_glob,
                'vp'        :vp_glob,
                'rosp'      :rosp_glob,
                'vortp'     :vortp_glob,
                'regime'    :regime_glob,
                'partsize'  :partsize_glob,
                'stokes'    :stokes_glob}

    if (ipart!=0 and igrowth ==0):
        return {'ug'        :ug_glob, 
                'vg'        :vg_glob,
                'rosg'      :rosg_glob,
                'pg'        :pg_glob,
                'vortg'     :vortg_glob,
                'divvg'     :divvg_glob,
                'up'        :up_glob,
                'vp'        :vp_glob,
                'rosp'      :rosp_glob,
                'vortp'     :vortp_glob}

    if (ipart==0 and igrowth ==0):
        return {'ug'        :ug_glob, 
                'vg'        :vg_glob,
                'rosg'      :rosg_glob,
                'pg'        :pg_glob,
                'vortg'     :vortg_glob,
                'divvg'     :divvg_glob}


def read_multiple_files(how_many_timesteps, arguments):
    def loopover(i,arguments):
        dictionary  = rf.read_single_timestep(i, arguments)

        new_dictionary['ug_%s' %i]          = dictionary['ug']
        new_dictionary['vg_%s' %i]          = dictionary['ug']
        new_dictionary['rosg_%s' %i]        = dictionary['ug']
        new_dictionary['pg_%s' %i]          = dictionary['ug']
        new_dictionary['vortg_%s' %i]       = dictionary['ug']
        new_dictionary['divvg_%s' %i]       = dictionary['ug']
        if ipart != 0:
            new_dictionary['up_%s' %i]          = dictionary['ug']
            new_dictionary['vp_%s' %i]          = dictionary['ug']
            new_dictionary['rosp_%s' %i]        = dictionary['ug']
            new_dictionary['vortp_%s' %i]       = dictionary['ug']
            if igrowth != 0:
                new_dictionary['regime_%s' %i]      = dictionary['ug']
                new_dictionary['partsize_%s' %i]    = dictionary['ug']
                new_dictionary['stokes_%s' %i]      = dictionary['ug']

        rpf.plot_singlelineplot(loop, arguments, data_stat['r'], data_stat['theta'], data['rosp'])



    new_dictionary = {}

    num_cores = multiprocessing.cpu_count()
    Parallel(n_jobs=num_cores)(delayed(loopover)(i, arguments) for i in range(1,4))


    for i in how_many_timesteps:
        dictionary = read_single_timestep(i, arguments)

        new_dictionary['ug_%s' %i]          = dictionary['ug']
        new_dictionary['vg_%s' %i]          = dictionary['ug']
        new_dictionary['rosg_%s' %i]        = dictionary['ug']
        new_dictionary['pg_%s' %i]          = dictionary['ug']
        new_dictionary['vortg_%s' %i]       = dictionary['ug']
        new_dictionary['divvg_%s' %i]       = dictionary['ug']
        if ipart != 0:
            new_dictionary['up_%s' %i]          = dictionary['ug']
            new_dictionary['vp_%s' %i]          = dictionary['ug']
            new_dictionary['rosp_%s' %i]        = dictionary['ug']
            new_dictionary['vortp_%s' %i]       = dictionary['ug']
            if igrowth != 0:
                new_dictionary['regime_%s' %i]      = dictionary['ug']
                new_dictionary['partsize_%s' %i]    = dictionary['ug']
                new_dictionary['stokes_%s' %i]      = dictionary['ug']
    
    return new_dictionary





