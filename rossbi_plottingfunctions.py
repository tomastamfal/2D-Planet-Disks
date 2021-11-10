from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib import rcParams
from pylab import *
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm
from matplotlib import colors
import cmocean 
import matplotlib.ticker


import rossbi_functions as rf

def cutforplot(matrix, arguments):
    nghost    	= arguments['nghost']
    new_matrix 	= matrix[nghost:-nghost, nghost:-nghost]
    return new_matrix

def calccellarea(r, plotting_arguments):
    # Cell area in cm^2
    AU_in_cm    = 1.49597870691*(10**(13))
    r           = r * plotting_arguments['R_ref'] 
    cell_area   = np.zeros(len(r))
    x_r         = np.zeros(len(r)+1)
    ratio       = (plotting_arguments['rout']/plotting_arguments['rin'])**(1.0/(plotting_arguments['nx_nodes']*plotting_arguments['xmax']))
    dtheta      = (2.0*np.pi)/(plotting_arguments['ny_nodes']* plotting_arguments['ymax'])*(180.0/np.pi)
    
    x_r[0] = r[0]/np.sqrt(ratio)
    for i in range(len(r)):
	x_r[i+1] = r[i]*np.sqrt(ratio)

    for i in range(len(r)):
        cell_area[i] = np.pi * (dtheta/360.0)* (x_r[i+1]**2.0 - r[i]**2.0) * (AU_in_cm**2.0)

    return cell_area



def plot_single2D_4times(outerloop, plotting_arguments, r, theta, data1, data2, data3, data4, plotrange=[]):
    rcParams.update({'font.size': plotting_arguments['fontsize']}) 
    bg_color    = plotting_arguments['bg_color']
    fg_color    = plotting_arguments['fg_color']
    label_color = plotting_arguments['label_color']
    colormap1   = plotting_arguments['colormap1']
    colormap2   = plotting_arguments['colormap2']
    colormap3   = plotting_arguments['colormap3']
    line_width  = plotting_arguments['line_width']
    name        = plotting_arguments['name']
    nghost 	= plotting_arguments['nghost']
    R_ref 	= plotting_arguments['R_ref']


    # define the colormap
    cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist[0] = (.5,.5,.5,1.0)
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    cmap = colors.ListedColormap(['black', 'red', 'aqua', 'yellow'])
    bounds = np.linspace(1,5,5)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)



    plot2D1    = cutforplot(data1, plotting_arguments)
    plot2D2    = cutforplot(data2, plotting_arguments)
    plot2D3    = cutforplot(data3, plotting_arguments)
    plot2D4    = cutforplot(data4, plotting_arguments)

    r_glob                       = cutforplot(r, plotting_arguments)*R_ref
    ytheta                       = np.linspace(0,2*np.pi, len(theta)-2*nghost) 
    radius_matrix , theta_matrix = np.meshgrid(r_glob[:,1],ytheta)
    
    if plotrange != []:
        vmin1 = plotrange[0]
        vmin2 = plotrange[2] 
        vmin3 = plotrange[4]
        
        vmax1 = plotrange[1] 
        vmax2 = plotrange[3]
        vmax3 = plotrange[5]


    else: 
        vmin1 = plot2D1.min()
        vmin2 = plot2D2.min()
        vmin3 = plot2D3.min()
        
        vmax1 = plot2D1.max()
        vmax2 = plot2D2.max()
        vmax3 = plot2D3.max()
    
    
    labelscbar1 = np.linspace(np.floor(np.log10(vmin1)), np.ceil(np.log10(vmax1)),num=1+(np.ceil(np.log10(vmax1))- np.floor(np.log10(vmin1))) )
    labelscbar2 = np.linspace(np.floor(np.log10(vmin2)), np.ceil(np.log10(vmax2)),num=1+(np.ceil(np.log10(vmax2))- np.floor(np.log10(vmin2))) )
    labelscbar3 = np.linspace(np.floor(np.log10(vmin3)), np.ceil(np.log10(vmax3)),num=1+(np.ceil(np.log10(vmax3))- np.floor(np.log10(vmin3))) )




    fig = plt.figure(figsize=(2*3.2677165354 , 2*2.3),dpi=600, facecolor=bg_color, edgecolor=fg_color)
    gs = GridSpec(2, 4, width_ratios=[1,1,1,1], height_ratios=[15,1])
    gs.update(left=0.01, right=0.99, hspace=0.0, wspace=0.15)

    ax1 = plt.subplot(gs[1,0], projection='polar', axisbg=bg_color)
    ax1.set_frame_on(False)
    ax1.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax1.spines.values():
            spine.set_color(fg_color)
    ima = ax1.pcolormesh(theta_matrix, radius_matrix,np.log10(plot2D1.transpose()), 
                    vmin=np.floor(np.log10(vmin1)), 
                    vmax=np.ceil(np.log10(vmax1)),cmap=colormap1)
    ax1.set_yticklabels([])
    ax1.set_theta_zero_location('N')
    thetaticks = []
    ax1.set_thetagrids(thetaticks, frac=1.3)   

    axes1 = plt.subplot(gs[0,0], axisbg=bg_color)
    cba = plt.colorbar(ima,cax=axes1, orientation="horizontal", ticks=labelscbar1)
    cba.ax.yaxis.set_ticks_position('left')
    cba.ax.yaxis.set_tick_params(color=fg_color, labelcolor=fg_color)  
    ax1.set_title(r"log$\left(\Sigma_d / \Sigma_{d,0}\right)$ ",color=fg_color)    
    cba.outline.set_linewidth(0.2)

    
    ax2 = plt.subplot(gs[1,1], projection='polar', axisbg=bg_color)
    ax2.set_frame_on(False)
    ax2.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax2.spines.values():
            spine.set_color(fg_color)
    ima = ax2.pcolormesh(theta_matrix, radius_matrix,np.log10(plot2D2.transpose()), 
                    vmin=np.floor(np.log10(vmin2)), 
                    vmax=np.ceil(np.log10(vmax2)),cmap=cmocean.cm.ice)
                    #vmax=np.ceil(np.log10(vmax2)),cmap=colormap2)
    ax2.set_yticklabels([])
    ax2.set_theta_zero_location('N')
    thetaticks = []
    ax2.set_thetagrids(thetaticks, frac=1.3)   

    axes2 = plt.subplot(gs[0,1], axisbg=bg_color)
    cba = plt.colorbar(ima,cax=axes2, orientation="horizontal", ticks=labelscbar2)
    cba.ax.yaxis.set_ticks_position('right')
    cba.ax.yaxis.set_tick_params(color=fg_color, labelcolor=fg_color)  
    ax2.set_title(r"log$\left(\Sigma_g/ \Sigma_{g,0} \right )$ ",color=fg_color)	    
    cba.outline.set_linewidth(0.2)

    ax3 = plt.subplot(gs[1,2], projection='polar', axisbg=bg_color)
    ax3.set_frame_on(False)
    ax3.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax3.spines.values():
            spine.set_color(fg_color)
    ima = ax3.pcolormesh(theta_matrix, radius_matrix,np.log10(plot2D3.transpose()), 
                    vmin=np.floor(np.log10(vmin3)), 
                    vmax=np.ceil(np.log10(vmax3)),cmap=colormap3)
    ax3.set_yticklabels([])
    ax3.set_theta_zero_location('N')
    thetaticks = []
    ax3.set_thetagrids(thetaticks, frac=1.3)   

    axes3 = plt.subplot(gs[0,2], axisbg=bg_color)
    cba = plt.colorbar(ima,cax=axes3, orientation="horizontal", ticks=labelscbar3)
    cba.ax.yaxis.set_ticks_position('left')
    cba.ax.yaxis.set_tick_params(color=fg_color, labelcolor=fg_color)  
    ax3.set_title(r"log$ \left (a_{\mathrm{max}} \right )$ [cm]",color=fg_color)    
    cba.outline.set_linewidth(0.2)
    
    print plot2D3.min()
    print plot2D3.max()
    
    ax4 = plt.subplot(gs[1,3], projection='polar', axisbg=bg_color)
    ax4.set_frame_on(False)
    ax4.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax4.yaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax4.spines.values():
            spine.set_color(fg_color)
    ima = ax4.pcolormesh(theta_matrix, radius_matrix,plot2D4.transpose(), 
                    vmin=0, 
                    vmax=5,cmap=cmap,norm=norm)
    ax4.set_yticklabels([])
    ax4.set_theta_zero_location('N')
    thetaticks = []
    ax4.set_thetagrids(thetaticks, frac=1.3)   

    axes4 = plt.subplot(gs[0,3], axisbg=bg_color)
    cba = plt.colorbar(ima,cax=axes4, orientation="horizontal")
    cba.ax.yaxis.set_ticks_position('right')
    labels = np.arange(0,5,1)
    loc = labels + .5
    cba.set_ticks(loc)
    cba.set_ticklabels(["frag","drift","df","a0"] )
    cba.ax.yaxis.set_tick_params(color=fg_color, labelcolor=fg_color)
    ax4.set_title(r"regimes",color=fg_color)	
    cba.outline.set_linewidth(0.2)

    plt.savefig('%04d_%s.png' %(outerloop, name), facecolor=bg_color, edgecolor=fg_color)
    plt.close()
    


def plot_single2D_2times(outerloop, plotting_arguments, r, theta, data1, data2, plotrange=[]):
    rcParams.update({'font.size': plotting_arguments['fontsize']}) 
    bg_color    = plotting_arguments['bg_color']
    fg_color    = plotting_arguments['fg_color']
    label_color = plotting_arguments['label_color']
    colormap1   = plotting_arguments['colormap1']
    colormap2   = plotting_arguments['colormap2']
    line_width  = plotting_arguments['line_width']
    name        = plotting_arguments['name']
    nghost 	= plotting_arguments['nghost']
    R_ref 	= plotting_arguments['R_ref']


    plot2D1    = cutforplot(data1, plotting_arguments)
    plot2D2    = cutforplot(data2, plotting_arguments)


    r_glob                       = cutforplot(r, plotting_arguments)*R_ref
    ytheta                       = np.linspace(0,2*np.pi, len(theta)-2*nghost) 
    radius_matrix , theta_matrix = np.meshgrid(r_glob[:,1],ytheta)

    if plotrange != []:
        vmin1 = plotrange[0]
        vmin2 = plotrange[2] 
        
        vmax1 = plotrange[1] 
        vmax2 = plotrange[3]


    else: 
        vmin1 = plot2D1.min()
        vmin2 = plot2D2.min()
        
        vmax1 = plot2D1.max()
        vmax2 = plot2D2.max()

    labelscbar1 = np.linspace(np.floor(np.log10(vmin1)), np.ceil(np.log10(vmax1)),num=1+(np.ceil(np.log10(vmax1))- np.floor(np.log10(vmin1))) )
    labelscbar2 = np.linspace(np.floor(np.log10(vmin2)), np.ceil(np.log10(vmax2)),num=1+(np.ceil(np.log10(vmax2))- np.floor(np.log10(vmin2))) )


    fig = plt.figure(figsize=(3.2677165354 , 2.3),dpi=600, facecolor=bg_color, edgecolor=fg_color)
    gs = GridSpec(2, 2, width_ratios=[1,1], height_ratios=[15,1])
    gs.update(left=0.01, right=0.99, hspace=0.0, wspace=0.15)

    ax1 = plt.subplot(gs[0,0], projection='polar', axisbg=bg_color)
    ax1.set_frame_on(False)
    ax1.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax1.spines.values():
            spine.set_color(fg_color)
    ima = ax1.pcolormesh(theta_matrix, radius_matrix,np.log10(plot2D1.transpose()), 
                    vmin=np.floor(np.log10(vmin1)), 
                    vmax=np.ceil(np.log10(vmax1)),cmap=colormap1)
    #ax1.set_rlabel_position(135)
    #rticks = [r_glob.min(), r_glob.max()]
    #rticks = []
    #ax1.set_rgrids(rticks, angle=300)
    ax1.set_yticklabels([])
    ax1.set_theta_zero_location('N')
    thetaticks = []
    ax1.set_thetagrids(thetaticks, frac=1.3)   

    axes1 = plt.subplot(gs[1,0], axisbg=bg_color)
    cba = plt.colorbar(ima,cax=axes1, orientation="horizontal", ticks=labelscbar1)
    cba.ax.yaxis.set_ticks_position('left')
    cba.set_ticks([-16, -11, -6, -1, 4 ])
    cba.ax.yaxis.set_tick_params(color=fg_color, labelcolor=fg_color, width=0.2)  
    ax1.set_title(r"log$ \left ( \Sigma_d / \Sigma_{d,0} \right)$ ",color=fg_color)    
    #cba.outline.set_visible(False)
    cba.outline.set_linewidth(0.2)

    
    ax2 = plt.subplot(gs[0,1], projection='polar', axisbg=bg_color)
    ax2.set_frame_on(False)
    #ax2.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax2.spines.values():
            spine.set_color(fg_color)
    ima = ax2.pcolormesh(theta_matrix, radius_matrix,np.log10(plot2D2.transpose()), 
                    vmin=np.floor(np.log10(vmin2)), 
                    vmax=np.ceil(np.log10(vmax2)),cmap=cmocean.cm.ice)
    #ax2.set_rlabel_position(135)
    #ax2.set_rgrids([r_glob.min(), r_glob.max()], angle=300)
    ax2.set_yticklabels([])
    ax2.set_theta_zero_location('N')

    thetaticks = []
    ax2.set_thetagrids(thetaticks, frac=1.3)   


    axes2 = plt.subplot(gs[1,1], axisbg=bg_color)
    cba = plt.colorbar(ima,cax=axes2, orientation="horizontal", ticks=labelscbar2)
    cba.ax.yaxis.set_ticks_position('right')
    cba.ax.yaxis.set_tick_params(color=fg_color, labelcolor=fg_color, width=0.2)  
    ax2.set_title(r"log$ \left ( \Sigma_g / \Sigma_{g,0}\right)$ ",color=fg_color)
    #cba.outline.set_visible(False)
    cba.outline.set_linewidth(0.2)
    
    plt.savefig('%04d_%s.png' %(outerloop, name), facecolor=bg_color, edgecolor=fg_color)
    plt.close()
    

def plot_singlelineplot(outerloop, plotting_arguments, r, data1, 
        plotting_arguments2=[], r2=[], data2=[], 
        plotting_arguments3=[], r3=[], data3=[], 
        plotting_arguments4=[], r4=[], data4=[], 
        plotting_arguments5=[], r5=[], data5=[], 
        planetpos=[],plotrange=[]):    
    rcParams.update({'font.size': plotting_arguments['fontsize']})    
    bg_color    = plotting_arguments['bg_color']
    fg_color    = plotting_arguments['fg_color']
    label_color = plotting_arguments['label_color']
    line_width  = plotting_arguments['line_width']
    name        = plotting_arguments['name']
    nghost 	= plotting_arguments['nghost']
    value_name  = plotting_arguments['value_name']
    value_label = plotting_arguments['value_label']
    R_ref 	= plotting_arguments['R_ref']
    ymax        = plotting_arguments['ymax']*plotting_arguments['ny_nodes']



    plot1    = cutforplot(data1, plotting_arguments)
    tmp1_line  = np.sum(plot1/(ymax), axis=1)

    if data2 != []:
        value_label2    = plotting_arguments2['value_label']
        ymax2           = plotting_arguments2['ymax']*plotting_arguments2['ny_nodes']
        plot2           = cutforplot(data2, plotting_arguments2)
        tmp2_line       = np.sum(plot2/(ymax2), axis=1)
        r_glob2         = cutforplot(r2, plotting_arguments2)
        r2              = r_glob2[:,1]*plotting_arguments2['R_ref']


    if data3 != []:
        value_label3    = plotting_arguments3['value_label']
        ymax3           = plotting_arguments3['ymax']*plotting_arguments3['ny_nodes']
        plot3           = cutforplot(data3, plotting_arguments3)
        tmp3_line       = np.sum(plot3/(ymax3), axis=1)
        r_glob3         = cutforplot(r3, plotting_arguments3)
        r3              = r_glob3[:,1]*plotting_arguments3['R_ref']


    if data4 != []:
        value_label4    = plotting_arguments4['value_label']
        ymax4           = plotting_arguments4['ymax']*plotting_arguments4['ny_nodes']
        plot4           = cutforplot(data4, plotting_arguments4)
        tmp4_line       = np.sum(plot4/(ymax4), axis=1)
        r_glob4         = cutforplot(r4, plotting_arguments4)
        r4              = r_glob4[:,1]*plotting_arguments4['R_ref']

    if data5 != []:
        value_label5    = plotting_arguments5['value_label']
        ymax5           = plotting_arguments5['ymax']*plotting_arguments5['ny_nodes']
        plot5           = cutforplot(data5, plotting_arguments5)
        tmp5_line       = np.sum(plot5/(ymax5), axis=1)
        r_glob5         = cutforplot(r5, plotting_arguments5)
        r5              = r_glob5[:,1]*plotting_arguments5['R_ref']


    r_glob                       = cutforplot(r, plotting_arguments)
    r                            = r_glob[:,1]*R_ref
 

    fig = plt.figure(figsize=(3.2677165354 , 3.3),dpi=600)
    gs = GridSpec(1, 1)
    gs.update(left=0.18, right=0.99, hspace=0.5, wspace=0.2, top=0.8, bottom=0.15)

    ax1 = plt.subplot(gs[0,0], axisbg=bg_color)
    ax1.set_frame_on(True)
    ax1.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax1.spines.values():
            spine.set_color(fg_color)
    ima = ax1.plot(r, tmp1_line,'k-', label='%s' %value_label)
    if data2 != []:
        ax1.plot(r2, tmp2_line,'r--', label='%s' %value_label2)
    if data3 != []:
        ax1.plot(r3, tmp3_line,'g-', label='%s' %value_label3)
    if data4 != []:
        ax1.plot(r4, tmp4_line,'b--', label='%s' %value_label4)
    if data5 != []:
        ax1.plot(r5, tmp5_line,'b-', label='%s' %value_label5)
    if planetpos != []:
        plt.axvline(planetpos, label= 'planet',color='r', linestyle='--', lw=2)
    
    if plotrange != []:
        plt.xlim([plotrange[0], plotrange[1]])
        plt.ylim([plotrange[2], plotrange[3]])

    # Put a legend below current axis
    ax1.legend(loc='upper center', bbox_to_anchor=(0.45, 1.3),
          fancybox=True, shadow=False, ncol = 3)
    
    plt.ylabel(str(value_name))
    plt.xlabel('Radius [AU]')

    plt.savefig('%04d_lineplot_%s.png' %(outerloop, name), facecolor=bg_color, edgecolor=fg_color)

def plot_singlelineplotzhu(outerloop, plotting_arguments, r, data1, 
        plotting_arguments2=[], r2=[], data2=[], 
        plotting_arguments3=[], r3=[], data3=[], 
        plotting_arguments4=[], r4=[], data4=[], 
        plotting_arguments5=[], r5=[], data5=[], 
        planetpos=[],plotrange=[]):    
    rcParams.update({'font.size': plotting_arguments['fontsize']}) 
    bg_color    = plotting_arguments['bg_color']
    fg_color    = plotting_arguments['fg_color']
    label_color = plotting_arguments['label_color']
    line_width  = plotting_arguments['line_width']
    name        = plotting_arguments['name']
    nghost 	= plotting_arguments['nghost']
    value_name  = plotting_arguments['value_name']
    value_label = plotting_arguments['value_label']
    R_ref 	= plotting_arguments['R_ref']
    ymax        = plotting_arguments['ymax']*plotting_arguments['ny_nodes']



    plot1    = cutforplot(data1, plotting_arguments)
    tmp1_line  = np.sum(plot1/(ymax), axis=1)
    value_label4    = plotting_arguments4['value_label']
    if data2 != []:
        value_label2    = plotting_arguments2['value_label']
        ymax2           = plotting_arguments2['ymax']*plotting_arguments2['ny_nodes']
        plot2           = cutforplot(data2, plotting_arguments2)
        tmp2_line       = np.sum(plot2/(ymax2), axis=1)
        r_glob2         = cutforplot(r2, plotting_arguments2)
        r2              = r_glob2[:,1]*plotting_arguments2['R_ref']


    if data3 != []:
        value_label3    = plotting_arguments3['value_label']
        ymax3           = plotting_arguments3['ymax']*plotting_arguments3['ny_nodes']
        plot3           = cutforplot(data3, plotting_arguments3)
        tmp3_line       = np.sum(plot3/(ymax3), axis=1)
        r_glob3         = cutforplot(r3, plotting_arguments3)
        r3              = r_glob3[:,1]*plotting_arguments3['R_ref']


    if data5 != []:
        value_label5    = plotting_arguments5['value_label']
        ymax5           = plotting_arguments5['ymax']*plotting_arguments5['ny_nodes']
        plot5           = cutforplot(data5, plotting_arguments5)
        tmp5_line       = np.sum(plot5/(ymax5), axis=1)
        r_glob5         = cutforplot(r5, plotting_arguments5)
        r5              = r_glob5[:,1]*plotting_arguments5['R_ref']


    r_glob                       = cutforplot(r, plotting_arguments)
    r                            = r_glob[:,1]*R_ref
 

    fig = plt.figure(figsize=(1.2*4.62677165354 , 1.8*3.3),dpi=600)
    gs = GridSpec(1, 1)
    gs.update(left=0.18, right=0.99, hspace=0.5, wspace=0.2, top=0.8, bottom=0.15)

    ax1 = plt.subplot(gs[0,0], axisbg=bg_color)
    ax1.set_frame_on(True)
    ax1.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax1.spines.values():
            spine.set_color(fg_color)
    ima = ax1.plot(r, tmp1_line,'g-', label='%s' %value_label)
    if data2 != []:
        ax1.plot(r2, tmp2_line,'r-', label='%s' %value_label2)
    if data3 != []:
        ax1.plot(r3, tmp3_line,'k-', label='%s' %value_label3)
    if data4 != []:
        ax1.plot(r4, data4 ,'b--', label='%s' %value_label4)
    if data5 != []:
        ax1.plot(r5, tmp5_line,'b-', label='%s' %value_label5)
    if planetpos != []:
        plt.axvline(planetpos, label= 'planet',color='r', linestyle='--', lw=2)
    
    if plotrange != []:
        plt.xlim([plotrange[0], plotrange[1]])
        plt.ylim([plotrange[2], plotrange[3]])

    # Put a legend below current axis
    #ax1.legend(loc='upper center', bbox_to_anchor=(0.385, 1.3),
    #      fancybox=True, shadow=False, ncol = 4)
    plt.legend()
    plt.ylabel(str(value_name))
    plt.xlabel('Radius [AU]')
    plt.axvspan(9.99,10.01, alpha=0.5, color='y')
    #plt.savefig('%04d_lineplot_%s.png' %(outerloop, name), facecolor=bg_color, edgecolor=fg_color)
    plt.savefig('%04d_lineplot_%s.pdf' %(outerloop, name))



def plot_singlehist(outerloop, plotting_arguments1, r1, data1_dust, data1_gas, data1_dust_norm=[], data1_gas_norm=[],
                    plotting_arguments2=[], r2=[], data2_dust=[], data2_gas=[], data2_dust_norm=[], data2_gas_norm=[],  
                    plotting_arguments3=[], r3=[], data3_dust=[], data3_gas=[], data3_dust_norm=[], data3_gas_norm=[], 
                    plotting_arguments4=[], r4=[], data4_dust=[], data4_gas=[], data4_dust_norm=[], data4_gas_norm=[], 
                    plotting_arguments5=[], r5=[], data5_dust=[], data5_gas=[], data5_dust_norm=[], data5_gas_norm=[], 
                    planetpos=[],plotrangedust=[], plotrangegas=[]):
    rcParams.update({'font.size': plotting_arguments1['fontsize']}) 
    g_in_Earthmass  = 1.6744252*(10**(-28))
    bg_color    = plotting_arguments1['bg_color']
    fg_color    = plotting_arguments1['fg_color']
    label_color = plotting_arguments1['label_color']
    name        = plotting_arguments1['name']
    nghost 	= plotting_arguments1['nghost']
    R_ref 	= plotting_arguments1['R_ref']
    ymax1       = plotting_arguments1['ymax']*plotting_arguments1['ny_nodes']
    ymax2       = plotting_arguments2['ymax']*plotting_arguments2['ny_nodes']
    ymax3       = plotting_arguments3['ymax']*plotting_arguments3['ny_nodes']
    ymax4       = plotting_arguments4['ymax']*plotting_arguments4['ny_nodes']
    ymax5       = plotting_arguments5['ymax']*plotting_arguments5['ny_nodes']

    color1  = 'black'
    color2  = 'red'
    color3  = 'blue'
    color4  = 'green'
    color5  = 'magenta'
    lw      = plotting_arguments1['line_width']


    r_glob1     = cutforplot(r1, plotting_arguments1)
    cell_area1  = calccellarea(r_glob1[:,1], plotting_arguments1)


    plot1_dust      = cutforplot(data1_dust, plotting_arguments1)
    plot1_gas       = cutforplot(data1_gas,  plotting_arguments1)
   
    tmpcell         = np.tile(cell_area1, (ymax1,1)).transpose()
    hist1_dust      = plot1_dust  *tmpcell *g_in_Earthmass
    hist1_gas       = plot1_gas   *tmpcell *g_in_Earthmass

    if data1_dust_norm != []:
        data1_dust_norm = cutforplot(data1_dust_norm,  plotting_arguments1)
        tmp1_dust       = np.sum(plot1_dust/(data1_dust_norm * ymax1), axis=1)
    else: 
      tmp1_dust       = np.sum(plot1_dust/(ymax1), axis=1)   
    
    if data1_gas_norm != []:
        data1_gas_norm = cutforplot(data1_gas_norm,  plotting_arguments1)
        tmp1_gas       = np.sum(plot1_gas/(data1_gas_norm * ymax1), axis=1)
    else: 
        tmp1_gas       = np.sum(plot1_gas/(ymax1), axis=1)   
    
    numberforbox = 1

    # FIND MIN AND MAX y axis
    plot_dust_min   = tmp1_dust.min()
    plot_dust_max   = tmp1_dust.max()
    plot_gas_min    = tmp1_gas.min()
    plot_gas_max    = tmp1_gas.max()
  
    hist_dust_min   = hist1_dust.min()
    hist_dust_max   = hist1_dust.max()
    hist_gas_min    = hist1_gas.min()
    hist_gas_max    = hist1_gas.max()


    if data2_dust != []:
        numberforbox += 1 
        r_glob2     = cutforplot(r2, plotting_arguments2)
        cell_area2  = calccellarea(r_glob2[:,1], plotting_arguments2)

        plot2_dust      = cutforplot(data2_dust, plotting_arguments2)
        plot2_gas       = cutforplot(data2_gas,  plotting_arguments2)
        
        tmpcell         = np.tile(cell_area2, (ymax2,1)).transpose()
        hist2_dust      = plot2_dust *tmpcell*g_in_Earthmass
        hist2_gas       = plot2_gas  *tmpcell*g_in_Earthmass

        if data2_dust_norm != []:
            data2_dust_norm = cutforplot(data2_dust_norm,  plotting_arguments2)
            tmp2_dust       = np.sum(plot2_dust/(data2_dust_norm * ymax2), axis=1)
        else: 
            tmp2_dust       = np.sum(plot2_dust/(ymax2), axis=1)
        
        if data2_gas_norm != []:
            data2_gas_norm = cutforplot(data2_gas_norm,  plotting_arguments2)
            tmp2_gas       = np.sum(plot2_gas/(data2_gas_norm * ymax2), axis=1)
        else: 
            tmp2_gas       = np.sum(plot2_gas/(ymax2), axis=1)
        
    
        # FIND MIN AND MAX
        plot_dust_min   = min(plot_dust_min, tmp2_dust.min())
        plot_dust_max   = max(plot_dust_max, tmp2_dust.max())
        plot_gas_min    = min(plot_gas_min,  tmp2_gas.min())
        plot_gas_max    = max(plot_gas_max,  tmp2_gas.max())
        
        hist_dust_min   = min(hist_dust_min, hist2_dust.min())
        hist_dust_max   = max(hist_dust_max, hist2_dust.max())
        hist_gas_min    = min(hist_gas_min,  hist2_gas.min())
        hist_gas_max    = max(hist_gas_max,  hist2_gas.max())



    if data3_dust != []:
        numberforbox += 1 

        r_glob3     = cutforplot(r3, plotting_arguments3)
        cell_area3  = calccellarea(r_glob3[:,1], plotting_arguments3)

        plot3_dust      = cutforplot(data3_dust, plotting_arguments3)
        plot3_gas       = cutforplot(data3_gas,  plotting_arguments3)
    
        tmpcell         = np.tile(cell_area3, (ymax3,1)).transpose()
        hist3_dust      = plot3_dust *tmpcell*g_in_Earthmass
        hist3_gas       = plot3_gas  *tmpcell*g_in_Earthmass

        if data3_dust_norm != []:
            data3_dust_norm = cutforplot(data3_dust_norm,  plotting_arguments3)
            tmp3_dust       = np.sum(plot3_dust/(data3_dust_norm * ymax3), axis=1)
        else: 
            tmp3_dust       = np.sum(plot3_dust/(ymax3), axis=1)
        
        if data3_gas_norm != []:
            data3_gas_norm = cutforplot(data3_gas_norm,  plotting_arguments3)
            tmp3_gas       = np.sum(plot3_gas/(data3_gas_norm * ymax3), axis=1)
        else: 
            tmp3_gas       = np.sum(plot3_gas/(ymax3), axis=1)
    
        # FIND MIN AND MAX
        plot_dust_min   = min(plot_dust_min, tmp3_dust.min())
        plot_dust_max   = max(plot_dust_max, tmp3_dust.max())
        plot_gas_min    = min(plot_gas_min,  tmp3_gas.min())
        plot_gas_max    = max(plot_gas_max,  tmp3_gas.max())
        
        hist_dust_min   = min(hist_dust_min, hist3_dust.min())
        hist_dust_max   = max(hist_dust_max, hist3_dust.max())
        hist_gas_min    = min(hist_gas_min,  hist3_gas.min())
        hist_gas_max    = max(hist_gas_max,  hist3_gas.max())

    if data4_dust != []:
        numberforbox += 1 

        r_glob4     = cutforplot(r4, plotting_arguments4)
        cell_area4  = calccellarea(r_glob4[:,1], plotting_arguments4)

        plot4_dust      = cutforplot(data4_dust, plotting_arguments4)
        plot4_gas       = cutforplot(data4_gas,  plotting_arguments4)

        tmpcell         = np.tile(cell_area4, (ymax4,1)).transpose()
        hist4_dust      = plot4_dust *tmpcell*g_in_Earthmass
        hist4_gas       = plot4_gas  *tmpcell*g_in_Earthmass

        if data4_dust_norm != []:
            data4_dust_norm = cutforplot(data4_dust_norm,  plotting_arguments4)
            tmp4_dust       = np.sum(plot4_dust/(data4_dust_norm * ymax4), axis=1)
        else: 
            tmp4_dust       = np.sum(plot4_dust/(ymax4), axis=1)
        
        if data4_gas_norm != []:
            data4_gas_norm = cutforplot(data4_gas_norm,  plotting_arguments4)
            tmp4_gas       = np.sum(plot4_gas/(data4_gas_norm * ymax4), axis=1)
        else: 
            tmp4_gas       = np.sum(plot4_gas/(ymax4), axis=1)
    
        # FIND MIN AND MAX
        plot_dust_min   = min(plot_dust_min, tmp4_dust.min())
        plot_dust_max   = max(plot_dust_max, tmp4_dust.max())
        plot_gas_min    = min(plot_gas_min,  tmp4_gas.min())
        plot_gas_max    = max(plot_gas_max,  tmp4_gas.max())
        
        hist_dust_min   = min(hist_dust_min, hist4_dust.min())
        hist_dust_max   = max(hist_dust_max, hist4_dust.max())
        hist_gas_min    = min(hist_gas_min,  hist4_gas.min())
        hist_gas_max    = max(hist_gas_max,  hist4_gas.max())


    if data5_dust != []:
        numberforbox += 1 

        r_glob5     = cutforplot(r5, plotting_arguments5)
        cell_area5  = calccellarea(r_glob5[:,1], plotting_arguments5)

        plot5_dust      = cutforplot(data5_dust, plotting_arguments5)
        plot5_gas       = cutforplot(data5_gas,  plotting_arguments5)


        tmpcell         = np.tile(cell_area5, (ymax5,1)).transpose()
        hist5_dust      = plot5_dust *tmpcell*g_in_Earthmass
        hist5_gas       = plot5_gas  *tmpcell*g_in_Earthmass

        if data5_dust_norm != []:
            data5_dust_norm = cutforplot(data5_dust_norm,  plotting_arguments5)
            tmp5_dust       = np.sum(plot5_dust/(data5_dust_norm * ymax5), axis=1)
        else: 
            tmp5_dust       = np.sum(plot5_dust/(ymax5), axis=1)
        
        if data5_gas_norm != []:
            data5_gas_norm = cutforplot(data5_gas_norm,  plotting_arguments5)
            tmp5_gas       = np.sum(plot5_gas/(data5_gas_norm * ymax5), axis=1)
        else: 
            tmp5_gas       = np.sum(plot5_gas/(ymax5), axis=1)
    
        # FIND MIN AND MAX
        plot_dust_min   = min(plot_dust_min, tmp5_dust.min())
        plot_dust_max   = max(plot_dust_max, tmp5_dust.max())
        plot_gas_min    = min(plot_gas_min,  tmp5_gas.min())
        plot_gas_max    = max(plot_gas_max,  tmp5_gas.max())
        
        hist_dust_min   = min(hist_dust_min, hist5_dust.min())
        hist_dust_max   = max(hist_dust_max, hist5_dust.max())
        hist_gas_min    = min(hist_gas_min,  hist5_gas.min())
        hist_gas_max    = max(hist_gas_max,  hist5_gas.max())

    #BINS 
    bingas  = np.logspace(floor(log10(hist_gas_min)),  ceil(log10(hist_gas_max)),  100)
    bindust = np.logspace(floor(log10(hist_dust_min)), ceil(log10(hist_dust_max)), 100)

    r1 = r_glob1[:,1]*R_ref
    #=======================================
    #HISTOGRAME SINGLE ---- NORMALIZED!

    clf()
    fig = plt.figure(figsize=(2*3.2677165354 , 5),dpi=600,facecolor=bg_color, edgecolor=fg_color)
    gs = GridSpec(2, 2, width_ratios=[1,1], height_ratios=[1,1])
    gs.update(left=0.1, right=0.86, hspace=0.15, wspace=0.1)

    ax1 = plt.subplot(gs[0,0], axisbg=bg_color)
    ax1.set_frame_on(True)
    ax1.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax1.yaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax1.spines.values():
	    spine.set_color(fg_color)

    ax1.hist(hist1_gas.flatten(), ls = '-',     lw = lw, bins = bingas, histtype='step', color=color1, label = plotting_arguments1['value_label'])
    if data2_gas != []:   
        ax1.hist(hist2_gas.flatten(), ls = 'dotted', lw = lw, bins = bingas, histtype='step', color=color2, label = plotting_arguments2['value_label'])
    if data3_gas != []:   
        ax1.hist(hist3_gas.flatten(), ls = 'dashed', lw = lw, bins = bingas, histtype='step', color=color3, label = plotting_arguments3['value_label'])
    if data4_gas != []:   
        ax1.hist(hist4_gas.flatten(), ls = '-',      lw = lw, bins = bingas, histtype='step', color=color4, label = plotting_arguments4['value_label'])
    if data5_gas != []:   
        ax1.hist(hist5_gas.flatten(), ls = 'dashed', lw = lw, bins = bingas, histtype='step', color=color5, label = plotting_arguments5['value_label'])


    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel("Number of cells", color = fg_color)
    if plotrangedust != []:
        ax1.set_ylim(plotrangedust)

    ax1.set_xlim([10**(floor(log10(hist_gas_min))), 10**(ceil(log10(hist_gas_max)))])



    ax2 = plt.subplot(gs[1,0], axisbg=bg_color)
    ax2.set_frame_on(True)
    ax2.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax2.yaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax2.spines.values():
	    spine.set_color(fg_color)
 
    ax2.hist(hist1_dust.flatten(), ls = '-',     lw = lw, bins = bindust, histtype='step', color=color1, label = plotting_arguments1['value_label'])
    if data2_dust != []:   
        ax2.hist(hist2_dust.flatten(), ls = 'dotted', lw = lw, bins = bindust, histtype='step', color=color2, label = plotting_arguments2['value_label'])
    if data3_dust != []:   
        ax2.hist(hist3_dust.flatten(), ls = 'dashed', lw = lw, bins = bindust, histtype='step', color=color3, label = plotting_arguments3['value_label'])
    if data4_dust != []:   
        ax2.hist(hist4_dust.flatten(), ls = '-',      lw = lw, bins = bindust, histtype='step', color=color4, label = plotting_arguments4['value_label'])
    if data5_dust != []:   
        ax2.hist(hist5_dust.flatten(), ls = 'dashed', lw = lw, bins = bindust, histtype='step', color=color5, label = plotting_arguments5['value_label'])   

    ax2.set_xscale('log')
    ax2.set_yscale('log')

    ax2.set_xlabel(r"$M /  M_{\oplus}$")
    ax2.set_ylabel("Number of cells", color = fg_color)

    if plotrangegas != []:
        ax2.set_ylim(plotrangegas)
    ax2.set_xlim([10**(floor(log10(hist_dust_min))), 10**(ceil(log10(hist_dust_max)))])
    

    ax3 = plt.subplot(gs[0,1], axisbg=bg_color)
    ax3.set_frame_on(True)
    ax3.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax3.yaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax3.spines.values():
	    spine.set_color(fg_color)
    ax3.plot(r1, tmp1_gas ,color=color1, linestyle='-',lw = lw, label= plotting_arguments1['value_label'])
    if data2_dust != []:
        r2 = r_glob2[:,1]*R_ref
        ax3.plot(r2, tmp2_gas, color=color2, linestyle=':', lw = lw, label= plotting_arguments2['value_label'])
    if data3_dust != []:
        r3 = r_glob3[:,1]*R_ref
        ax3.plot(r3, tmp3_gas, color=color3, linestyle='--',lw = lw, label= plotting_arguments3['value_label'])
    if data4_dust != []:
        r4 = r_glob4[:,1]*R_ref
        ax3.plot(r4, tmp4_gas, color=color4, linestyle='-.',lw = lw, label= plotting_arguments4['value_label'])
    if data5_dust != []:
        r5 = r_glob5[:,1]*R_ref
        ax3.plot(r5, tmp5_gas, color=color5, linestyle='-', lw = lw, label= plotting_arguments5['value_label'])
        
    ax3.set_xscale('log')
    ax3.set_yscale('log')

    if data1_gas_norm != []:
        ax3.set_ylabel(r'$\sum_{i=1}^N \Sigma_{g,i}/\Sigma_0 /N $')
    else:
        ax3.set_ylabel(r'$\sum_{i= 1}^N \Sigma_{g,i} /N$')

    ax3.set_ylim([10**(floor(log10(plot_gas_min -0.1* plot_gas_min))), 10**(ceil(log10(plot_gas_max+0.1*plot_gas_max)))])
    ax3.set_xlim([np.floor(r1.min()), np.around(r1.max())])
    ax3.set_xticks([np.floor(r1.min()), R_ref, np.around(r1.max())])    
    ax3.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax3.xaxis.set_ticklabels(['5','','','','','10','','','','','','','','','','','','','','','','','','','','30'],minor=True)
    ax3.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax3.yaxis.set_label_position("right")
    ax3.yaxis.tick_right()

    ax4 = plt.subplot(gs[1,1],axisbg=bg_color)
    ax4.set_frame_on(True)
    for spine in ax4.spines.values():
	    spine.set_color(fg_color)

    ax4.plot(r1, tmp1_dust, color=color1, linestyle='-', lw = lw, label= plotting_arguments1['value_label']) 
    
    if data2_dust != []:   
        ax4.plot(r2, tmp2_dust, color=color2, linestyle=':',  lw = lw, label= plotting_arguments2['value_label'])
    if data3_dust != []:   
        ax4.plot(r3, tmp3_dust, color=color3, linestyle='--', lw = lw, label= plotting_arguments3['value_label'])
    if data4_dust != []:   
        ax4.plot(r4, tmp4_dust, color=color4, linestyle='-.', lw = lw, label= plotting_arguments4['value_label'])
    if data5_dust != []:   
        ax4.plot(r5, tmp5_dust, color=color5, linestyle='-',  lw = lw, label= plotting_arguments5['value_label'])
    
    if data1_gas_norm != []:
        ax4.set_ylabel(r'$\sum_{i=1}^N \Sigma_{d,i}/\Sigma_0 /N $')
    else:
        ax4.set_ylabel(r'$\sum_{i=1}^N \Sigma_{d,i}/N $')
    ax4.set_ylim([10**(floor(log10(plot_dust_min))), 10**(ceil(log10(plot_dust_max)))])

    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_xlabel('Radius [AU]', color = fg_color)
    ax4.set_xlim([np.floor(r1.min()), np.around(r1.max())])
    ax4.set_xticks([np.floor(r1.min()), R_ref, np.around(r1.max())])  
    ax4.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax4.xaxis.set_ticklabels(['5','','','','','10','','','','','','','','','','','','','','','','','','','','30'],minor=True)
    ax4.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax4.yaxis.set_label_position("right")
    ax4.yaxis.tick_right()


    # Put a legend below current axis
    ax4.legend(loc='upper center', bbox_to_anchor=(0, 2.4),
          fancybox=True, shadow=False, ncol=numberforbox)
    plt.savefig('%04d_compare_%s.png' %(outerloop, name), facecolor=bg_color, edgecolor=fg_color)
    plt.close()



def plot_compare_2D(outerloop, plotting_arguments, r, theta, data1, data2, data3, data4, data5, plotting_arguments2, r2, theta2, data6, data7, data8, plotrange=[]):
    
    bg_color    = plotting_arguments['bg_color']
    fg_color    = plotting_arguments['fg_color']
    label_color = plotting_arguments['label_color']
    colormap1   = plotting_arguments['colormap1']
    colormap4   = plotting_arguments['colormap2']
    colormap3   = plotting_arguments['colormap3']
    line_width  = plotting_arguments['line_width']
    name        = plotting_arguments['name']
    nghost1 	= plotting_arguments['nghost']
    nghost2 	= plotting_arguments2['nghost']
    R_ref1 	= plotting_arguments['R_ref']
    R_ref2 	= plotting_arguments2['R_ref']

    fixedsize   = plotting_arguments2['a']

    rcParams.update({'font.size': plotting_arguments['fontsize']}) 

    plot2D1    = cutforplot(data1, plotting_arguments)
    plot2D2    = cutforplot(data2, plotting_arguments)
    plot2D3    = cutforplot(data3, plotting_arguments)
    plot2D4    = cutforplot(data4, plotting_arguments)
    plot2D5    = cutforplot(data5, plotting_arguments)
    plot2D6    = cutforplot(data6, plotting_arguments2)
    plot2D7    = cutforplot(data7, plotting_arguments2)
    plot2D8    = cutforplot(data8, plotting_arguments2)
    plot2D9    = plot2D8*0.0 + fixedsize

    r_glob1                        = cutforplot(r, plotting_arguments)*R_ref1
    ytheta1                        = np.linspace(0,2*np.pi, len(theta)-2*nghost1) 
    radius_matrix1 , theta_matrix1 = np.meshgrid(r_glob1[:,1],ytheta1)

    r_glob2                        = cutforplot(r2, plotting_arguments2)*R_ref2
    ytheta2                        = np.linspace(0,2*np.pi, len(theta2)-2*nghost2) 
    radius_matrix2 , theta_matrix2 = np.meshgrid(r_glob2[:,1],ytheta2)

    if plotrange != []:
        vmin1 = plotrange[0]
        vmin2 = plotrange[2] 
        vmin3 = plotrange[4] 
        vmin4 = min(plot2D4.min(), fixedsize)
        vmin4 = max(vmin4, 1e-5)        
        
        vmax1 = plotrange[1] 
        vmax2 = plotrange[3]
        vmax3 = plotrange[5] 
        vmax4 = max(plot2D4.max(),fixedsize)


    else: 
        vmin1 = min(plot2D1.min(), plot2D6.min())
        vmin2 = min(plot2D2.min(), plot2D7.min())
        vmin3 = min(plot2D3.min(), plot2D8.min())
        vmin4 = min(plot2D4.min(), fixedsize)
        vmin4 = max(vmin4, 1e-5)
       
        vmax1 = max(plot2D1.max(), plot2D6.max())
        vmax2 = max(plot2D2.max(), plot2D7.max())
        vmax3 = max(plot2D3.max(), plot2D8.max())
        vmax4 = max(plot2D4.max(), fixedsize)

    # define the colormap
    cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist[0] = (.5,.5,.5,1.0)
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    cmap = colors.ListedColormap(['black', 'red', 'aqua', 'yellow'])
    bounds = np.linspace(1,5,5)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    print vmin4
    
#    labelscbar1 = np.linspace(np.floor(np.log10(vmin1)), np.ceil(np.log10(vmax1)),num=1+(np.ceil(np.log10(vmax1))- np.floor(np.log10(vmin1))) )
#    labelscbar2 = np.linspace(np.floor(np.log10(vmin2)), np.ceil(np.log10(vmax2)),num=1+(np.ceil(np.log10(vmax2))- np.floor(np.log10(vmin2))) )
#    labelscbar3 = np.linspace(np.floor(np.log10(vmin3)), np.ceil(np.log10(vmax3)),num=1+(np.ceil(np.log10(vmax3))- np.floor(np.log10(vmin3))) )
    labelscbar1 = np.linspace(round(np.floor(np.log10(vmin1)),2), round(np.ceil(np.log10(vmax1)),2), num=3)
    labelscbar2 = np.linspace(round(np.floor(np.log10(vmin2)),1), round(np.ceil(np.log10(vmax2)),1), num=3)   
    labelscbar4 = np.linspace(round(np.floor(np.log10(vmin4)),1), round(np.ceil(np.log10(vmax4)),1), num=3)      
    labelscbar3 = np.linspace(round(np.floor(np.log10(vmin3)),1), round(np.ceil(np.log10(vmax3)),1), num=3)    
#    labelscbar4 = np.linspace(np.floor(np.log10(vmin4)), np.ceil(np.log10(vmax4)),num=1+(np.ceil(np.log10(vmax4))- np.floor(np.log10(vmin4)) ))
    labelscbar1 = [round(np.log10(vmin1),1), 0, round(np.log10(vmax1),1)]
    print labelscbar1
    fig = plt.figure(figsize=(4*3.2677165354 , 3*2.3), facecolor=bg_color, edgecolor=fg_color, dpi=150)
    gs = GridSpec(3, 5, width_ratios=[1,1,1,1,1], height_ratios=[15,15,1])
    gs.update(left=0.01, right=0.99, hspace=0.0, wspace=0.15)


# =============
# BOTTOM LINE
# =============


    ax1 = plt.subplot(gs[1,0], projection='polar', axisbg=bg_color)
    ax1.set_frame_on(False)
    ax1.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax1.spines.values():
            spine.set_color(fg_color)
    ima = ax1.pcolormesh(theta_matrix1, radius_matrix1,np.log10(plot2D1.transpose()), 
                    #vmin=np.floor(np.log10(vmin1)),
                    vmin=round(np.log10(vmin1),1),
                    #vmax=np.ceil(np.log10(vmax1)),
                    vmax=round(np.log10(vmax1),1),
                    cmap=colormap1)
    ax1.set_yticklabels([])
    ax1.set_theta_zero_location('N')
    thetaticks = []
    ax1.set_thetagrids(thetaticks, frac=1.3)   

    
    ax2 = plt.subplot(gs[1,1], projection='polar', axisbg=bg_color)
    ax2.set_frame_on(False)
    ax2.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax2.spines.values():
            spine.set_color(fg_color)
    ima = ax2.pcolormesh(theta_matrix1, radius_matrix1,np.log10(plot2D2.transpose()), 
                    #vmin=np.floor(np.log10(vmin2)),
                    vmin=round(np.log10(vmin2),1),
                    #vmax=np.ceil(np.log10(vmax2)),
                    vmax=round(np.log10(vmax2),1),
                    cmap=cmocean.cm.ice)
                    #vmax=np.ceil(np.log10(vmax2)),cmap=colormap2)
    ax2.set_yticklabels([])
    ax2.set_theta_zero_location('N')
    thetaticks = []
    ax2.set_thetagrids(thetaticks, frac=1.3)   



    ax3 = plt.subplot(gs[1,2], projection='polar', axisbg=bg_color)
    ax3.set_frame_on(False)
    ax3.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax3.spines.values():
            spine.set_color(fg_color)
    ima = ax3.pcolormesh(theta_matrix1, radius_matrix1,np.log10(plot2D3.transpose()), 
                    #vmin=np.floor(np.log10(vmin2)),
                    vmin=round(np.floor(np.log10(vmin3)),1),
                    #vmax=np.ceil(np.log10(vmax2)),
                    vmax=round(np.ceil(np.log10(vmax3)),1),
                    cmap=colormap3)
    ax3.set_yticklabels([])
    ax3.set_theta_zero_location('N')
    thetaticks = []
    ax3.set_thetagrids(thetaticks, frac=1.3)   



    ax4 = plt.subplot(gs[1,3], projection='polar', axisbg=bg_color)
    ax4.set_frame_on(False)
    ax4.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax4.spines.values():
            spine.set_color(fg_color)
    ima = ax4.pcolormesh(theta_matrix1, radius_matrix1, np.log10(plot2D4.transpose()), 
                    #vmin=np.floor(np.log10(vmin2)),
                    vmin=round(np.floor(np.log10(vmin4)),1),
                    #vmax=np.ceil(np.log10(vmax2)),
                    vmax=round(np.ceil(np.log10(vmax4)),1),
                    cmap=colormap4)
                    #vmax=np.ceil(np.log10(vmax2)),cmap=colormap2)
    ax4.set_yticklabels([])
    ax4.set_theta_zero_location('N')
    thetaticks = []
    ax4.set_thetagrids(thetaticks, frac=1.3)  




    ax5 = plt.subplot(gs[1,4], projection='polar', axisbg=bg_color)
    ax5.set_frame_on(False)
    ax5.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax5.yaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax5.spines.values():
            spine.set_color(fg_color)
    ima = ax5.pcolormesh(theta_matrix1, radius_matrix1,plot2D5.transpose(), 
                    vmin=0, 
                    vmax=5,cmap=cmap,norm=norm)
    ax5.set_yticklabels([])
    ax5.set_theta_zero_location('N')
    thetaticks = []
    ax5.set_thetagrids(thetaticks, frac=1.3)   


    axes5 = plt.subplot(gs[2,4], axisbg=bg_color)
    cba = plt.colorbar(ima,cax=axes5, orientation="horizontal")
    cba.ax.yaxis.set_ticks_position('right')
    labels = np.arange(0,5,1)
    loc = labels + .5
    cba.set_ticks(loc)
    cba.set_ticklabels(["frag","drift","df","a0"] )
    cba.ax.yaxis.set_tick_params(color=fg_color, labelcolor=fg_color)
    ax5.set_title(r"regimes",color=fg_color)	
    cba.outline.set_linewidth(0.2)


# =============
#TOP LINE
# =============


    ax6 = plt.subplot(gs[0,0], projection='polar', axisbg=bg_color)
    ax6.set_frame_on(False)
    ax6.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax6.spines.values():
            spine.set_color(fg_color)
    ima = ax6.pcolormesh(theta_matrix2, radius_matrix2,np.log10(plot2D6.transpose()), 
                    #vmin=np.floor(np.log10(vmin1)),
                    vmin=round(np.log10(vmin1),1),
                    #vmax=np.ceil(np.log10(vmax1)),
                    vmax=round(np.log10(vmax1),1),
                    cmap=colormap1)
    ax6.set_yticklabels([])
    ax6.set_theta_zero_location('N')
    thetaticks = []
    ax6.set_thetagrids(thetaticks, frac=1.3)   

    axes6 = plt.subplot(gs[2,0], axisbg=bg_color)
    cba = plt.colorbar(ima,cax=axes6, orientation="horizontal", ticks=labelscbar1)
    cba.ax.yaxis.set_ticks_position('left')
    cba.ax.yaxis.set_tick_params(color=fg_color, labelcolor=fg_color)  
    ax6.set_title(r"log$\left(\Sigma_d / \Sigma_{d,0}\right)$ ",color=fg_color)    
    cba.outline.set_linewidth(0.2)


    ax7 = plt.subplot(gs[0,1], projection='polar', axisbg=bg_color)
    ax7.set_frame_on(False)
    ax7.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax7.spines.values():
            spine.set_color(fg_color)
    ima = ax7.pcolormesh(theta_matrix2, radius_matrix2,np.log10(plot2D7.transpose()), 
                    #vmin=np.floor(np.log10(vmin2)),
                    vmin=round(np.log10(vmin2),1),
                    #vmax=np.ceil(np.log10(vmax2)),
                    vmax=round(np.log10(vmax2),1),
                    cmap=cmocean.cm.ice)
                    #vmax=np.ceil(np.log10(vmax2)),cmap=colormap2)
    ax7.set_yticklabels([])
    ax7.set_theta_zero_location('N')
    thetaticks = []
    ax7.set_thetagrids(thetaticks, frac=1.3)  

    axes7 = plt.subplot(gs[2,1], axisbg=bg_color)
    cba = plt.colorbar(ima,cax=axes7, orientation="horizontal", ticks=labelscbar2)
    cba.ax.yaxis.set_ticks_position('right')
    cba.ax.yaxis.set_tick_params(color=fg_color, labelcolor=fg_color)  
    ax7.set_title(r"log$\left(\Sigma_g/ \Sigma_{g,0} \right )$ ",color=fg_color)	    
    cba.outline.set_linewidth(0.2)
    


    ax8 = plt.subplot(gs[0,2], projection='polar', axisbg=bg_color)
    ax8.set_frame_on(False)
    ax8.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax8.spines.values():
            spine.set_color(fg_color)
    ima = ax8.pcolormesh(theta_matrix2, radius_matrix2, np.log10(plot2D8.transpose()), 
                    #vmin=np.floor(np.log10(vmin2)),
                    vmin=round(np.ceil(np.log10(vmin3)),1),
                    #vmax=np.ceil(np.log10(vmax2)),
                    vmax=round(np.floor(np.log10(vmax3)),1),
                    cmap=colormap3)
                    #vmax=np.ceil(np.log10(vmax2)),cmap=colormap2)
    ax8.set_yticklabels([])
    ax8.set_theta_zero_location('N')
    thetaticks = []
    ax8.set_thetagrids(thetaticks, frac=1.3)   

    axes8 = plt.subplot(gs[2,2], axisbg=bg_color)
    cba = plt.colorbar(ima,cax=axes8, orientation="horizontal", ticks=labelscbar3)
    cba.ax.yaxis.set_ticks_position('right')
    cba.ax.yaxis.set_tick_params(color=fg_color, labelcolor=fg_color)  
    ax8.set_title(r" log$(\mathrm{St})$",color=fg_color)	    
    cba.outline.set_linewidth(0.2)


    ax9 = plt.subplot(gs[0,3], projection='polar', axisbg=bg_color)
    ax9.set_frame_on(False)
    ax9.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax9.spines.values():
            spine.set_color(fg_color)
    ima = ax9.pcolormesh(theta_matrix2, radius_matrix2,np.log10(plot2D9.transpose()), 
                    #vmin=np.floor(np.log10(vmin2)),
                    vmin=round(np.floor(np.log10(vmin4)),1),
                    #vmax=np.ceil(np.log10(vmax2)),
                    vmax=round(np.ceil(np.log10(vmax4)),1),
                    cmap=colormap4)
    ax9.set_yticklabels([])
    ax9.set_theta_zero_location('N')
    thetaticks = []
    ax9.set_thetagrids(thetaticks, frac=1.3)  

    axes9 = plt.subplot(gs[2,3], axisbg=bg_color)
    cba = plt.colorbar(ima,cax=axes9, orientation="horizontal", ticks=labelscbar4)
    cba.ax.yaxis.set_ticks_position('left')
    cba.ax.yaxis.set_tick_params(color=fg_color, labelcolor=fg_color)  
    ax9.set_title(r"log(a$_{max}$) [cm]",color=fg_color)	
    cba.outline.set_linewidth(0.2)



    plt.savefig('%04d_%s_new.png' %(outerloop, name), facecolor=bg_color, edgecolor=fg_color, dpi=150)
    #plt.savefig('%04d_%s.eps' %(outerloop, name), format='eps')
    #plt.savefig('%04d_%s.pdf' %(outerloop, name) )



    plt.close()


def single_analysis(outerloop, r1, data1_dust, data1_gas, plotting_arguments1, data1_dust_norm=[], data1_gas_norm=[], plotrange=[]):
    rcParams.update({'font.size': plotting_arguments['fontsize']})     
    bg_color    = plotting_arguments1['bg_color']
    fg_color    = plotting_arguments1['fg_color']
    label_color = plotting_arguments1['label_color']
    name        = plotting_arguments1['name']
    nghost 	= plotting_arguments1['nghost']
    R_ref 	= plotting_arguments1['R_ref']
    ymax1       = plotting_arguments1['ymax']*plotting_arguments1['ny_nodes']
    
    g_in_Earthmass  = 1.6744252*(10**(-28))
    
    color1  = 'black'
    color2  = 'red'
    color3  = 'blue'
    color4  = 'green'
    color5  = 'magenta'
    lw      = plotting_arguments1['line_width']


    r_glob1     = cutforplot(r1, plotting_arguments1)
    cell_area1  = calccellarea(r_glob1[:,1], plotting_arguments1)


    plot1_dust      = cutforplot(data1_dust, plotting_arguments1)
    plot1_gas       = cutforplot(data1_gas,  plotting_arguments1)

    hist1_dust      = plot1_dust *cell_area1*g_in_Earthmass
    hist1_gas       = plot1_gas  *cell_area1*g_in_Earthmass

    if data1_dust_norm != []:
        data1_dust_norm = cutforplot(data1_dust_norm, plotting_arguments1)      
        tmp1_dust       = np.sum(plot1_dust/(data1_dust_norm * ymax1), axis=1)
    else: 
      tmp1_dust         = np.sum(plot1_dust/(ymax1), axis=1)   
    
    if data1_gas_norm != []:
        data1_gas_norm  = cutforplot(data1_gas_norm, plotting_arguments1) 
        tmp1_gas        = np.sum(plot1_gas/(data1_gas_norm * ymax1), axis=1)
    else: 
        tmp1_gas        = np.sum(plot1_gas/(ymax1), axis=1)   


    # FIND MIN AND MAX y axis
    plot_dust_min   = tmp1_dust.min()
    plot_dust_max   = tmp1_dust.max()
    plot_gas_min    = tmp1_gas.min()
    plot_gas_max    = tmp1_gas.max()
  
    hist_dust_min   = hist1_dust.min()
    hist_dust_max   = hist1_dust.max()
    hist_gas_min    = hist1_gas.min()
    hist_gas_max    = hist1_gas.max()

    #BINS 
    bingas  = np.logspace(floor(log10(hist_gas_min)),  ceil(log10(hist_gas_max)),  100)
    bindust = np.logspace(floor(log10(hist_dust_min)), ceil(log10(hist_dust_max)), 100)

    r1 = r_glob1[:,1]*R_ref
    #=======================================
    #HISTOGRAME SINGLE ---- NORMALIZED!

    clf()
    fig = plt.figure(figsize=(30,25),facecolor=bg_color, edgecolor=fg_color)
    gs = GridSpec(2, 2, width_ratios=[1,1], height_ratios=[1,1])
    gs.update(left=0.15, right=0.9, hspace=0.15, wspace=0.5)

    ax1 = plt.subplot(gs[1,0], axisbg=bg_color)
    ax1.set_frame_on(True)
    ax1.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax1.yaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax1.spines.values():
	    spine.set_color(fg_color)
    ax1.hist(hist1_dust.flatten(), lw = lw, bins = bindust, histtype='step', color=color2, label = 'dust')
    ax1.set_xscale('log')
    ax1.set_ylabel("Number of cells", color = fg_color)
    if plotrange != []:
        ax1.set_ylim(plotrangedust)
    plt.legend(loc=1) 



    ax2 = plt.subplot(gs[1,1], axisbg=bg_color)
    ax2.set_frame_on(True)
    ax2.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax2.yaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax2.spines.values():
	    spine.set_color(fg_color)
    
    ax2.hist(hist1_gas.flatten(), lw = lw, bins = bingas, histtype='step', color=color1, label = 'gas')
    ax2.set_xscale('log')
    ax2.set_xlabel(r"$M /  M_{\bigoplus} $"+"\n\n", color = fg_color)    
    ax2.set_ylabel("Number of cells", color = fg_color)
    if plotrange != []:
        ax2.set_ylim(plotrangegas)
    plt.legend(loc=1) 

    ax3 = plt.subplot(gs[0,:], axisbg=bg_color)
    ax3.set_frame_on(True)
    ax3.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax3.yaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax3.spines.values():
	    spine.set_color(fg_color)
    ax3.plot(r1, tmp1_gas,  color=color1, linestyle='--', lw = lw, label= 'gas')
    ax3.plot(r1, tmp1_dust, color=color2, linestyle=':', lw = lw, label= 'dust') 
    ax3.set_yscale('log')
    if data1_gas_norm != []:
        ax3.set_ylabel(r'$\sum_{i = 1}^{N} \Sigma_i / \Sigma_0/N$', color = fg_color)

    else:
        ax3.set_ylabel(r'$\sum_{i = 1}^{N} \Sigma_i / N$', color = fg_color)

    plotmin = min(plot_gas_min, plot_dust_min)
    plotmax = max(plot_gas_max, plot_dust_max) 
    ax3.set_ylim([10**(floor(log10(plotmin -0.1* plotmin))), 10**(ceil(log10(plotmax+0.1*plotmax)))])
    ax3.set_xlim([np.floor(r1.min()), np.around(r1.max())])
    ax3.set_xticks([np.floor(r1.min()), R_ref, np.around(r1.max())])    
    ax3.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.legend(loc=4) 
    plt.savefig('%04d_analysis_%s.png' %(outerloop, name), facecolor=bg_color, edgecolor=fg_color)





def plot_flux(outerloop,r, theta, planetdata, data, x_vel, y_vel, plotting_arguments, Hill_arguments):
    rcParams.update({'font.size': plotting_arguments['fontsize']})     
    bg_color    = plotting_arguments['bg_color']
    fg_color    = plotting_arguments['fg_color']
    label_color = plotting_arguments['label_color']
    colormap1   = plotting_arguments['colormap1']
    colormap2   = plotting_arguments['colormap2']
    line_width  = plotting_arguments['line_width']
    name        = plotting_arguments['name']
    nghost 	= plotting_arguments['nghost']
    R_ref 	= plotting_arguments['R_ref']
    r_in 	= plotting_arguments['r_in']
    r_out 	= plotting_arguments['r_out']
    nx          = plotting_arguments['xmax']*plotting_arguments['nx_nodes']
    ny          = plotting_arguments['ymax']*plotting_arguments['ny_nodes']

    #TIME CONVERS.
    AU   = 1.5e13                           # IN [cm]     
    year = 3.154e7                          # IN [s]
    v_ref  = R_ref**(-0.5) * 2.97396e6      # IN [cm/s]8 s_ref  = 2*np.pi*R_ref*AU           # IN [cm]
    s_ref  = 2*np.pi*R_ref*AU               # IN [cm]
     
    time = s_ref/v_ref                      # IN [s]
    time = time/year                        # IN [yr]

    multiple_plot_radius = 3.0



    x_planet    = planetdata[0::10,1]*R_ref
    y_planet    = planetdata[0::10,2]
    time_planet = planetdata[0::10,0]
    vx_planet   = planetdata[0::10,3]*v_ref
    vy_planet   = planetdata[0::10,4]*v_ref



    plot2D      = cutforplot(data, plotting_arguments)
    plot2Dx     = cutforplot(x_vel, plotting_arguments)
    plot2Dy     = cutforplot(y_vel, plotting_arguments)

    r_glob      = cutforplot(r, plotting_arguments)
    theta_glob  = theta[nghost:-nghost]

    ratio    = (r_out/r_in)**(1.0/(nx))
    x_c      = r_glob[:,1]*R_ref
    x        = np.zeros(nx+1)

    planet_radius   = Hill_arguments['planet_radius']
    multiple        = Hill_arguments['multiple']
    m               = Hill_arguments['m']
    M               = Hill_arguments['M']

    planet_radius   = 0.00046732617 #AU
    multiple        = 5000
    r               = multiple  * planet_radius

    # HILL RADIUS:
    m  = 1.898e27 #kg 
    M  = 1.989e30 #kg

    if (Hill_arguments['Hillradius'] == 'True'):
        print '=============================='
        print 'Using Hill radius:'
        r = R_ref*(m/(3.0*M))**(1.0/3.0)
        print r
    else:
        r = multiple  * planet_radius
        print '=============================='
        print 'Using multiple m and radius is:'
        print r


    if (r < 4*(R_ref*(r_out-r_in))/nx ):
        print '=============================='
        print 'Change of the radius!!!!!!'
        print '=============================='
        print 'Old r = %s' %r
        r = 4*(R_ref*(r_out-r_in))/nx
        print 'New r = %s' %r



    #GRID IN POLAR COORDINATES!!!!!
    dx      = np.zeros(nx)
    x[0]    = x_c[1]/np.sqrt(ratio)  
    for i in range(nx):
        x[i+1]    = x_c[i] * np.sqrt(ratio)
        dx[i]     = x[i+1]-x[i]




    dy   = 2*np.pi/(ny)
    y_c  = theta_glob
    print np.shape(y_c)
    y    = np.zeros(ny+1)
    for j in range(ny):
        y[j+1] = y_c[j] + 0.5*dy

    y[0] = y_c[0] - 0.5*dy



    # TRANSFORMATION TO CARTESIAN!!!!
    # CENTERS!
    x_cartesian = np.zeros(nx*ny)
    y_cartesian = np.zeros(nx*ny)
    circle      = np.zeros([nx,ny])
    velocity    = np.zeros([nx,ny])
    print np.shape(circle)


    k = 0
    for i in range(nx):
        for j in range(ny):
            x_cartesian[k]  = x_c[i]* np.cos(y_c[j])
            y_cartesian[k]  = x_c[i]* np.sin(y_c[j])
            tmp_r           = np.sqrt((x_planet[int(outerloop)]-x_cartesian[k])**2 + (y_planet[int(outerloop)]-y_cartesian[k])**2)
            if (tmp_r < r ):
                velocity[i,j]   = np.sqrt((vx_planet[int(outerloop)]-plot2Dx[i,j])**2 + (vy_planet[int(outerloop)]-plot2Dy[i,j])**2)*v_ref
                tmp             = 2*np.pi* tmp_r *plot2D[i,j]
                circle[i,j]     = tmp *velocity[i,j]*3.154e+7*time/(5.9736e27)
            else:
                circle[i,j] = np.nan

            k += 1

    r_m, theta = np.meshgrid(x_c, y_c)

    X = r_m* np.cos(theta)
    Y = r_m* np.sin(theta)
    origin = 'lower'

    lev_exp = np.arange(np.floor(-2),np.ceil(2))
    levs = np.power(10, lev_exp)
    levs = [1e-1, 2e-1, 3e-1,4e-1,5e-1,6e-1,7e-1,8e-1,9e-1,1e0, 2e0, 3e0,4e0, 5e0]
    
    lev_exp2 = np.arange(np.floor(-10),np.ceil(-2))
    levs2 = np.power(10, lev_exp2)
    levs2 = [9e-6,1e-5,1.1e-5,1.2e-5,1.3e-5,1.4e-5,1.5e-5,1.6e-5,1.7e-5,1.8e-5,1.9e-5, 2e-5,2.1e-5,2.2e-5,2.3e-5,2.4e-5,2.5e-5,2.6e-5,2.7e-5,2.8e-5,2.9e-5, 3e-5]

    plt.figure(figsize=(15,8))
    #PLOT NEW CIRCLE
    my1_cmap = plt.cm.inferno
    my1_cmap.set_bad('white', 1.)
    C=plt.pcolormesh(X,Y, plot2D.transpose(), cmap=plt.cm.inferno_r)
    #C =plt.contourf(X,Y, plot2D.transpose(), levs,  cmap=my1_cmap,origin=origin)
    #C =plt.contourf(X,Y, plot2D.transpose(), cmap=my1_cmap,origin=origin, norm=matplotlib.colors.LogNorm()) 

    my2_cmap = plt.cm.ocean_r
    my2_cmap.set_bad('white', 1.)
    masked_array =np.ma.masked_where(circle == np.nan, circle) 
    #CS = plt.contourf(X,Y, masked_array.transpose(), levs2, cmap=my2_cmap, origin=origin)
    CS = plt.pcolormesh(X,Y,masked_array.transpose(),  cmap=my2_cmap)
    
    
    #my3_cmap = plt.cm.bone_r
    #my3_cmap.set_bad('white', 1.)
    #CS2 = plt.contour(X,Y,masked_array.transpose(), levels=CS.levels, cmap=my3_cmap,origin=origin, norm=colors.LogNorm())
    #CS2 = plt.contour(X,Y,masked_array.transpose(), levels=CS.levels, cmap=my3_cmap,origin=origin, norm=colors.LogNorm())

    #PLOT PLANET
    plt.scatter(x_planet[int(outerloop)], y_planet[int(outerloop)], color = 'red')
    
    cbar = plt.colorbar(CS)
    #cbar.add_lines(CS2)
    cbar.ax.set_ylabel(r'[M$_\oplus$/orbit]')

    cbar2 = plt.colorbar(C, orientation = 'horizontal')
    cbar2.ax.set_xlabel(r'Dust surface density [g/cm$^2$]')

    plt.xlabel('[AU]')
    plt.ylabel('[AU]')
    plt.xlim([x_planet[int(outerloop)] - multiple_plot_radius*r, x_planet[int(outerloop)]+multiple_plot_radius*r])
    plt.ylim([y_planet[int(outerloop)] - multiple_plot_radius*r, y_planet[int(outerloop)]+multiple_plot_radius*r])
    plt.savefig('%05d_FLUX.png' %outerloop)



def partsize_analysis(outerloop, plotting_arguments, r, theta, dust, gas, pressure,cs,time, planetpos=[],plotrange=[]): 
    rcParams.update({'font.size': plotting_arguments['fontsize']}) 
    # define the colormap
    cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist[0] = (.5,.5,.5,1.0)
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    cmap = colors.ListedColormap(['black', 'red', 'aqua', 'yellow'])
    bounds = np.linspace(1,5,5)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)


    bg_color    = plotting_arguments['bg_color']
    fg_color    = plotting_arguments['fg_color']
    label_color = plotting_arguments['label_color']
    line_width  = plotting_arguments['line_width']
    name        = plotting_arguments['name']
    nghost 	= plotting_arguments['nghost']
    R_ref 	= plotting_arguments['R_ref']
    sig_ref 	= plotting_arguments['sig_ref']
    v_ref 	= plotting_arguments['v_ref']
    ymax        = plotting_arguments['ymax']*plotting_arguments['ny_nodes']
    xmax        = plotting_arguments['xmax']*plotting_arguments['nx_nodes']
    dust_gas    = plotting_arguments['dust_gas']
    a0          = plotting_arguments['a0']
    rho_sol     = plotting_arguments['rho_sol']

    rho_p      = cutforplot(dust, plotting_arguments)
    rho_g      = cutforplot(gas,  plotting_arguments)
    p_g        = pressure
    cs_g       = cutforplot(cs, plotting_arguments)
    
    rloop      = r[:,1]
    r_glob     = cutforplot(r, plotting_arguments)
    r          = r_glob[:,1]*R_ref
    
    ytheta                       = np.linspace(0,2*np.pi, len(theta)-2*nghost) 
    radius_matrix , theta_matrix = np.meshgrid(r,ytheta)

    x_ref = R_ref*1.5e13 

    smallv     = 1.e-28            # ridiculously small value
    lalpha     = 1e-3              # turbulence strength parameter
    vf         = 1000              # fragmentation threshold [cm/s]

    afrag  = np.zeros((xmax, ymax))
    adrift = np.zeros((xmax, ymax))
    adf    = np.zeros((xmax, ymax))
    aini   = np.zeros((xmax, ymax))
   
    gammadust = np.zeros((xmax, ymax))
    test_matrix = np.zeros((xmax,ymax))
    regime = test_matrix
    omega_k = rloop[nghost:-nghost]**(-1.5)



    testungmax = -10000
    testungmin = 1e111
    for j in range(xmax):
        for i in range(ymax):
            
            gammadust[i,j] = abs(rloop[i+nghost] * (p_g[i+1+nghost,j+nghost]-p_g[i+nghost,j+nghost]) / (p_g[i+nghost,j+nghost] * (rloop[i+1+nghost]-rloop[i+nghost])))

            afrag[i,j]     =  (0.37 * 2.0 * (rho_g[i,j] *sig_ref) * vf**2.0) / (3.0 * np.pi * rho_sol * lalpha *(cs_g[i,j]  * v_ref)**2.0)
            adrift[i,j]    =  0.2 * (1.1 * (rho_p[i,j] *sig_ref) * (omega_k[i]*v_ref/x_ref * rloop[i]*x_ref)**2.) / (np.pi * rho_sol * (cs_g[i,j]  * v_ref)**2.0 * abs(gammadust[i,j]))
            adf[i,j]       =  0.2 * (2.0 * (rho_g[i,j] *sig_ref) * vf * (omega_k[i]*v_ref/x_ref * rloop[i]*x_ref)) / (np.pi * rho_sol * abs(gammadust[i,j]) * (cs_g[i,j]*v_ref)**2.0 * 0.5)
            aini[i,j]      =  a0 * np.exp(0.8 * time *(x_ref/v_ref)* omega_k[i]*(v_ref/x_ref)* dust_gas)
           
           

            test_matrix[i,j] = min(afrag[i,j], adrift[i,j], adf[i,j])
    
    	    tmp_rad                = test_matrix[i,j]
            test_matrix[i,j]       = min(aini[i,j] , test_matrix[i,j] )
		
            if ((afrag[i,j]   < adrift[i,j])   and (afrag[i,j]   < adf[i,j] )):
                regime[i,j]  = 1  
            if ((adrift[i,j]   < afrag[i,j])   and (adrift[i,j] < adf[i,j]  )):
                regime[i,j]  = 2  
            if ((adf[i,j]   < afrag[i,j])      and (adf[i,j]    < adrift[i,j]  )):
                regime[i,j]  = 3
            if (aini[i,j]   < tmp_rad):
                regime[i,j]  = 4 
   
            print test_matrix[i,j]
            if (testungmax < test_matrix[i,j]):
                testungmax = test_matrix[i,j]
            if (testungmin > test_matrix[i,j]):
                testungmin = test_matrix[i,j]

        print testungmax
        print testungmin
        print '----------'
    print '----------'
    print testungmax
    print testungmin
    print '----------'


    
    tmp1_line       = np.sum(afrag, axis=1)/ymax
    tmp2_line       = np.sum(adrift, axis=1)/ymax
    tmp3_line       = np.sum(adf, axis=1)/ymax
    tmp4_line       = np.sum(aini, axis=1)/ymax


    stokes1         = np.sum((np.pi *rho_sol*afrag/(2.*rho_g*sig_ref)) /ymax,axis=1) 
    stokes2         = np.sum((np.pi *rho_sol*adrift/(2.*rho_g*sig_ref))/ymax,axis=1)
    stokes3         = np.sum((np.pi *rho_sol*adf/(2.*rho_g*sig_ref))/ymax, axis=1)
    stokes4         = np.sum((np.pi *rho_sol*aini/(2.*rho_g*sig_ref))/ymax,axis=1)   
    
    
    fig = plt.figure(figsize=(3.2677165354 , 3.3),dpi=600)
    gs = GridSpec(1, 1)
    gs.update(left=0.18, right=0.99, hspace=0.5, wspace=0.2, top=0.8, bottom=0.15)

    ax1 = plt.subplot(gs[0,0], axisbg=bg_color)
    ax1.set_frame_on(True)
    ax1.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax1.spines.values():
            spine.set_color(fg_color)
   
    ima = ax1.loglog(r, tmp1_line,'k-', label=r'a$_{frag}$')
    ax1.loglog(r, tmp2_line,'r--', label=r'a$_{drift}$')
    ax1.loglog(r, tmp3_line,'g-', label=r'a$_{df}$')
    ax1.loglog(r, tmp4_line,'b--', label=r'a$_{0}$')
    
    if plotrange != []:
        plt.xlim([plotrange[0], plotrange[1]])
        plt.ylim([plotrange[2], plotrange[3]])
    
    ax1.set_ylim([1e-4,1e2])

    # Put a legend below current axis
    ax1.legend(loc='upper center', bbox_to_anchor=(0.45, 1.3),
          fancybox=True, shadow=False, ncol = 3)
    
    plt.ylabel('particle size')
    plt.xlabel('Radius [AU]')

    plt.savefig('%04d_partsize_%s.png' %(outerloop, name), facecolor=bg_color, edgecolor=fg_color)

    plt.close()


    fig = plt.figure(figsize=(3.2677165354 , 3.3),dpi=600)
    gs = GridSpec(1, 1)
    gs.update(left=0.18, right=0.99, hspace=0.5, wspace=0.2, top=0.8, bottom=0.15)

    ax1 = plt.subplot(gs[0,0], axisbg=bg_color)
    ax1.set_frame_on(True)
    ax1.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax1.spines.values():
            spine.set_color(fg_color)
   
    ima = ax1.loglog(r, stokes1,'k-', label=r'a$_{frag}$')
    ax1.loglog(r, stokes2,'r--', label=r'a$_{drift}$')
    ax1.loglog(r, stokes3,'g-', label=r'a$_{df}$')
    ax1.loglog(r, stokes4,'b--', label=r'a$_{0}$')
    
    if plotrange != []:
        plt.xlim([plotrange[0], plotrange[1]])
        plt.ylim([plotrange[2], plotrange[3]])
    
    ax1.set_ylim([1e-5,1])

    # Put a legend below current axis
    ax1.legend(loc='upper center', bbox_to_anchor=(0.45, 1.3),
          fancybox=True, shadow=False, ncol = 3)
    
    plt.ylabel('Stokes number')
    plt.xlabel('Radius [AU]')

    plt.savefig('%04d_partsize_stokes_%s.png' %(outerloop, name), facecolor=bg_color, edgecolor=fg_color)

    plt.close()


    fig = plt.figure(figsize=(30,27),facecolor=bg_color, edgecolor=fg_color)
    gs = GridSpec(2,1 , width_ratios=[10], height_ratios=[10,1])
    gs.update(left=0.1, right=0.9, hspace=0.0, wspace=0.0)


    ax4 = plt.subplot(gs[0,0], projection='polar', axisbg=bg_color)
    ax4.set_frame_on(False)
    ax4.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax4.yaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax4.spines.values():
            spine.set_color(fg_color)
    ima = ax4.pcolormesh(theta_matrix, radius_matrix, regime.transpose(), 
                    vmin=0, 
                    vmax=5,cmap=cmap,norm=norm)
    ax4.set_yticklabels([])
    ax4.set_theta_zero_location('N')
    thetaticks = []
    ax4.set_thetagrids(thetaticks, frac=1.3)   

    axes4 = plt.subplot(gs[1,0], axisbg=bg_color)
    cba = plt.colorbar(ima,cax=axes4, orientation="horizontal")
    cba.ax.yaxis.set_ticks_position('right')
    labels = np.arange(0,5,1)
    loc = labels + .5
    cba.set_ticks(loc)
    cba.set_ticklabels(["frag","drift","df","a0"] )
    cba.ax.yaxis.set_tick_params(color=fg_color, labelcolor=fg_color)
    ax4.set_title(r"regimes",color=fg_color)	
    cba.outline.set_linewidth(0.2)

    
    plt.savefig('Test.png')



def plot_paperhist(outerloop, plotting_arguments1, r1, data1_dust, data1_gas,   data1_partsize, data1_stokes, data1_dust_norm=[], data1_gas_norm=[],
                    plotting_arguments2=[], r2=[], data2_dust=[], data2_gas=[], data2_partsize=[], data2_stokes=[], data2_dust_norm=[], data2_gas_norm=[],  
                    plotting_arguments3=[], r3=[], data3_dust=[], data3_gas=[], data3_partsize=[], data3_stokes=[], data3_dust_norm=[], data3_gas_norm=[], 
                    plotting_arguments4=[], r4=[], data4_dust=[], data4_gas=[], data4_partsize=[], data4_stokes=[], data4_dust_norm=[], data4_gas_norm=[], 
                    plotting_arguments5=[], r5=[], data5_dust=[], data5_gas=[], data5_partsize=[], data5_stokes=[], data5_dust_norm=[], data5_gas_norm=[], 
                    planetpos=[],plotrangedust=[], plotrangegas=[]):
    rcParams.update({'font.size': plotting_arguments1['fontsize']}) 
    bg_color    = plotting_arguments1['bg_color']
    fg_color    = plotting_arguments1['fg_color']
    label_color = plotting_arguments1['label_color']
    name        = plotting_arguments1['name']
    nghost 	= plotting_arguments1['nghost']
    R_ref 	= plotting_arguments1['R_ref']
    amin 	= min(plotting_arguments1['a'],plotting_arguments2['a'],plotting_arguments3['a'], plotting_arguments4['a'],plotting_arguments5['a'])/10.0
    ymax1       = plotting_arguments1['ymax']*plotting_arguments1['ny_nodes']
    ymax2       = plotting_arguments2['ymax']*plotting_arguments2['ny_nodes']
    ymax3       = plotting_arguments3['ymax']*plotting_arguments3['ny_nodes']
    ymax4       = plotting_arguments4['ymax']*plotting_arguments4['ny_nodes']
    ymax5       = plotting_arguments5['ymax']*plotting_arguments5['ny_nodes']

    color1  = 'black'
    color2  = 'red'
    color3  = 'blue'
    color4  = 'green'
    color5  = 'magenta'
    lw      = plotting_arguments1['line_width']


    r_glob1     = cutforplot(r1, plotting_arguments1)


    plot1_dust      = cutforplot(data1_dust, plotting_arguments1)
    plot1_gas       = cutforplot(data1_gas,  plotting_arguments1)
    plot1_partsize  = cutforplot(data1_partsize,  plotting_arguments1)
    plot1_stokes    = cutforplot(data1_stokes,  plotting_arguments1)


    if data1_dust_norm != []:
        data1_dust_norm = cutforplot(data1_dust_norm,  plotting_arguments1)
        tmp1_dust       = np.sum(plot1_dust/(data1_dust_norm * ymax1), axis=1)
    else: 
      tmp1_dust       = np.sum(plot1_dust/(ymax1), axis=1)   
    
    if data1_gas_norm != []:
        data1_gas_norm = cutforplot(data1_gas_norm,  plotting_arguments1)
        tmp1_gas       = np.sum(plot1_gas/(data1_gas_norm * ymax1), axis=1)
    else: 
        tmp1_gas       = np.sum(plot1_gas/(ymax1), axis=1)   
    
    tmp1_partsize       = np.sum(plot1_partsize/(ymax1), axis=1) 
    tmp1_stokes         = np.sum(plot1_stokes/(ymax1), axis=1) 
    numberforbox = 1

    # FIND MIN AND MAX y axis
    plot_dust_min   = tmp1_dust.min()
    plot_dust_max   = tmp1_dust.max()
    plot_gas_min    = tmp1_gas.min()
    plot_gas_max    = tmp1_gas.max()
    
    
    plot_stokes_min    = tmp1_stokes.min()
    plot_stokes_max    = tmp1_stokes.max()
    plot_partsize_min  = tmp1_partsize.min()
    plot_partsize_max  = tmp1_partsize.max()


    if data2_dust != []:
        numberforbox += 1 
        r_glob2     = cutforplot(r2, plotting_arguments2)

        plot2_dust      = cutforplot(data2_dust, plotting_arguments2)
        plot2_gas       = cutforplot(data2_gas,  plotting_arguments2)
        plot2_partsize  = cutforplot(data2_partsize,  plotting_arguments2)
        plot2_stokes    = cutforplot(data2_stokes,  plotting_arguments2)


        if data2_dust_norm != []:
            data2_dust_norm = cutforplot(data2_dust_norm,  plotting_arguments2)
            tmp2_dust       = np.sum(plot2_dust/(data2_dust_norm * ymax2), axis=1)
        else: 
            tmp2_dust       = np.sum(plot2_dust/(ymax2), axis=1)
        
        if data2_gas_norm != []:
            data2_gas_norm = cutforplot(data2_gas_norm,  plotting_arguments2)
            tmp2_gas       = np.sum(plot2_gas/(data2_gas_norm * ymax2), axis=1)
        else: 
            tmp2_gas       = np.sum(plot2_gas/(ymax2), axis=1)
        
        tmp2_partsize       = np.sum(plot2_partsize/(ymax2), axis=1) 
        tmp2_stokes         = np.sum(plot2_stokes/(ymax2), axis=1) 
    
        # FIND MIN AND MAX
        plot_dust_min   = min(plot_dust_min, tmp2_dust.min())
        plot_dust_max   = max(plot_dust_max, tmp2_dust.max())
        plot_gas_min    = min(plot_gas_min,  tmp2_gas.min())
        plot_gas_max    = max(plot_gas_max,  tmp2_gas.max())
        
        plot_stokes_min   = min(plot_stokes_min, tmp2_stokes.min())
        plot_stokes_max   = max(plot_stokes_max, tmp2_stokes.max())
        plot_partsize_min = min(plot_partsize_min,  tmp2_partsize.min())
        plot_partsize_max = max(plot_partsize_max,  tmp2_partsize.max())


    if data3_dust != []:
        numberforbox += 1 

        r_glob3     = cutforplot(r3, plotting_arguments3)

        plot3_dust      = cutforplot(data3_dust, plotting_arguments3)
        plot3_gas       = cutforplot(data3_gas,  plotting_arguments3)
        plot3_partsize  = cutforplot(data3_partsize,  plotting_arguments3)
        plot3_stokes    = cutforplot(data3_stokes,  plotting_arguments3)

        if data3_dust_norm != []:
            data3_dust_norm = cutforplot(data3_dust_norm,  plotting_arguments3)
            tmp3_dust       = np.sum(plot3_dust/(data3_dust_norm * ymax3), axis=1)
        else: 
            tmp3_dust       = np.sum(plot3_dust/(ymax3), axis=1)
        
        if data3_gas_norm != []:
            data3_gas_norm = cutforplot(data3_gas_norm,  plotting_arguments3)
            tmp3_gas       = np.sum(plot3_gas/(data3_gas_norm * ymax3), axis=1)
        else: 
            tmp3_gas       = np.sum(plot3_gas/(ymax3), axis=1)

        tmp3_partsize       = np.sum(plot3_partsize/(ymax3), axis=1) 
        tmp3_stokes         = np.sum(plot3_stokes/(ymax3), axis=1) 

        # FIND MIN AND MAX
        plot_dust_min   = min(plot_dust_min, tmp3_dust.min())
        plot_dust_max   = max(plot_dust_max, tmp3_dust.max())
        plot_gas_min    = min(plot_gas_min,  tmp3_gas.min())
        plot_gas_max    = max(plot_gas_max,  tmp3_gas.max())

        plot_stokes_min   = min(plot_stokes_min, tmp3_stokes.min())
        plot_stokes_max   = max(plot_stokes_max, tmp3_stokes.max())
        plot_partsize_min = min(plot_partsize_min,  tmp3_partsize.min())
        plot_partsize_max = max(plot_partsize_max,  tmp3_partsize.max())

        
    if data4_dust != []:
        numberforbox += 1 

        r_glob4         = cutforplot(r4, plotting_arguments4)
        plot4_dust      = cutforplot(data4_dust, plotting_arguments4)
        plot4_gas       = cutforplot(data4_gas,  plotting_arguments4)
        plot4_partsize  = cutforplot(data4_partsize,  plotting_arguments4)
        plot4_stokes    = cutforplot(data4_stokes,  plotting_arguments4)

        if data4_dust_norm != []:
            data4_dust_norm = cutforplot(data4_dust_norm,  plotting_arguments4)
            tmp4_dust       = np.sum(plot4_dust/(data4_dust_norm * ymax4), axis=1)
        else: 
            tmp4_dust       = np.sum(plot4_dust/(ymax4), axis=1)
        
        if data4_gas_norm != []:
            data4_gas_norm = cutforplot(data4_gas_norm,  plotting_arguments4)
            tmp4_gas       = np.sum(plot4_gas/(data4_gas_norm * ymax4), axis=1)
        else: 
            tmp4_gas       = np.sum(plot4_gas/(ymax4), axis=1)
        
        tmp4_partsize       = np.sum(plot4_partsize/(ymax4), axis=1) 
        tmp4_stokes         = np.sum(plot4_stokes/(ymax4), axis=1) 

        # FIND MIN AND MAX
        plot_dust_min   = min(plot_dust_min, tmp4_dust.min())
        plot_dust_max   = max(plot_dust_max, tmp4_dust.max())
        plot_gas_min    = min(plot_gas_min,  tmp4_gas.min())
        plot_gas_max    = max(plot_gas_max,  tmp4_gas.max())
        
        plot_stokes_min   = min(plot_stokes_min, tmp4_stokes.min())
        plot_stokes_max   = max(plot_stokes_max, tmp4_stokes.max())
        plot_partsize_min = min(plot_partsize_min,  tmp4_partsize.min())
        plot_partsize_max = max(plot_partsize_max,  tmp4_partsize.max())


    if data5_dust != []:
        numberforbox += 1 

        r_glob5     = cutforplot(r5, plotting_arguments5)

        plot5_dust      = cutforplot(data5_dust, plotting_arguments5)
        plot5_gas       = cutforplot(data5_gas,  plotting_arguments5)
        plot5_partsize  = cutforplot(data5_partsize,  plotting_arguments5)
        plot5_stokes    = cutforplot(data2_stokes,  plotting_arguments5)


        if data5_dust_norm != []:
            data5_dust_norm = cutforplot(data5_dust_norm,  plotting_arguments5)
            tmp5_dust       = np.sum(plot5_dust/(data5_dust_norm * ymax5), axis=1)
        else: 
            tmp5_dust       = np.sum(plot5_dust/(ymax5), axis=1)
        
        if data5_gas_norm != []:
            data5_gas_norm = cutforplot(data5_gas_norm,  plotting_arguments5)
            tmp5_gas       = np.sum(plot5_gas/(data5_gas_norm * ymax5), axis=1)
        else: 
            tmp5_gas       = np.sum(plot5_gas/(ymax5), axis=1)
        
        tmp5_partsize       = np.sum(plot5_partsize/(ymax5), axis=1) 
        tmp5_stokes         = np.sum(plot5_stokes/(ymax5), axis=1) 

        # FIND MIN AND MAX
        plot_dust_min   = min(plot_dust_min, tmp5_dust.min())
        plot_dust_max   = max(plot_dust_max, tmp5_dust.max())
        plot_gas_min    = min(plot_gas_min,  tmp5_gas.min())
        plot_gas_max    = max(plot_gas_max,  tmp5_gas.max())
        
        plot_stokes_min   = min(plot_stokes_min, tmp5_stokes.min())
        plot_stokes_max   = max(plot_stokes_max, tmp5_stokes.max())
        plot_partsize_min = min(plot_partsize_min,  tmp5_partsize.min())
        plot_partsize_max = max(plot_partsize_max,  tmp5_partsize.max())

    
    #CUT THE PLOT
    plot_dust_min      = max(plot_dust_min, 1e-6)
    plot_partsize_max  = max(1e1, plot_partsize_max)

    r1 = r_glob1[:,1]*R_ref

    #=======================================
    #HISTOGRAME SINGLE ---- NORMALIZED!

    clf()
    fig = plt.figure(figsize=(2*3.2677165354 , 5),dpi=150,facecolor=bg_color, edgecolor=fg_color)
    gs = GridSpec(2, 2, width_ratios=[1,1], height_ratios=[1,1])
    gs.update(left=0.15, right=0.86, hspace=0.15, wspace=0.1)

    ax3 = plt.subplot(gs[1,0], axisbg=bg_color)
    ax3.set_frame_on(True)
    ax3.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax3.yaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax3.spines.values():
	    spine.set_color(fg_color)

    ax3.plot(r1, tmp1_partsize ,color=color1, linestyle='-',lw = lw, label= plotting_arguments1['value_label'])
    if data2_gas != []:
        r2 = r_glob2[:,1]*R_ref
        ax3.plot(r2, tmp2_partsize, ls = 'dotted', lw = lw, color=color2, label = plotting_arguments2['value_label'])
    if data3_gas != []:   
        r3 = r_glob3[:,1]*R_ref
        ax3.plot(r3, tmp3_partsize, ls = 'dashed', lw = lw, color=color3, label = plotting_arguments3['value_label'])
    if data4_gas != []:   
        r4 = r_glob4[:,1]*R_ref
        ax3.plot(r4, tmp4_partsize, ls = '-',      lw = lw, color=color4, label = plotting_arguments4['value_label'])
    if data5_gas != []:   
        r5 = r_glob5[:,1]*R_ref
        ax3.plot(r5, tmp5_partsize, ls = 'dashed', lw = lw, color=color5, label = plotting_arguments5['value_label'])


    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlabel('Radius [AU]', color = fg_color)
    ax3.set_ylabel(r"$\sum_{i=1}^N $ a$_{max,i}/N$ [cm]", color = fg_color)
    if plotrangedust != []:
        ax3.set_ylim(plotrangedust)

    ax3.set_ylim([min(plot_partsize_min,amin), max(plot_partsize_max,amin)])
    ax3.set_xlim([np.floor(r1.min()), np.around(r1.max())])
    ax3.set_xticks([np.floor(r1.min()), R_ref, np.around(r1.max())])    
    ax3.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax3.xaxis.set_ticklabels(['5','','','','','10','','','','','','','','','','','','','','','','','','','','30'],minor=True)
    ax3.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax3.yaxis.set_label_position("left")
    ax3.yaxis.tick_left()    



    ax4 = plt.subplot(gs[1,1], axisbg=bg_color)
    ax4.set_frame_on(True)
    ax4.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax4.yaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax4.spines.values():
	    spine.set_color(fg_color)
 
    ax4.plot(r1, tmp1_stokes ,color=color1, linestyle='-',lw = lw, label= plotting_arguments1['value_label'])
    if data2_gas != []:   
        r2 = r_glob2[:,1]*R_ref
        ax4.plot(r2, tmp2_stokes, ls = 'dotted', lw = lw, color=color2, label = plotting_arguments2['value_label'])
    if data3_gas != []:   
        r3 = r_glob3[:,1]*R_ref
        ax4.plot(r3, tmp3_stokes, ls = 'dashed', lw = lw, color=color3, label = plotting_arguments3['value_label'])
    if data4_gas != []:   
        r4 = r_glob4[:,1]*R_ref
        ax4.plot(r4, tmp4_stokes, ls = '-',      lw = lw, color=color4, label = plotting_arguments4['value_label'])
    if data5_gas != []:   
        r5 = r_glob5[:,1]*R_ref
        ax4.plot(r5, tmp5_stokes, ls = 'dashed', lw = lw, color=color5, label = plotting_arguments5['value_label'])


    ax4.set_xscale('log')
    ax4.set_yscale('log')

    ax4.set_xlabel('Radius [AU]', color = fg_color)
    ax4.set_ylabel(r"$\sum_{i=1}^N $ $\mathrm{St}_{i}/N$", color = fg_color)

    if plotrangegas != []:
        ax4.set_ylim(plotrangegas)
    #ax4.set_xlim([10**(floor(log10(plot_stokes_min))), 10**(ceil(log10(plot_stokes_max)))])
    #ax4.set_ylim([10**(floor(log10(plot_gas_min -0.1* plot_gas_min))), 10**(ceil(log10(plot_gas_max+0.1*plot_gas_max)))])
    ax4.set_xlim([np.floor(r2.min()), np.around(r2.max())])
    ax4.set_xticks([np.floor(r2.min()), R_ref, np.around(r2.max())])    
    ax4.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax4.xaxis.set_ticklabels(['5','','','','','10','','','','','','','','','','','','','','','','','','','','30'],minor=True)
    ax4.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax4.yaxis.set_label_position("right")
    ax4.yaxis.tick_right()    

    ax1 = plt.subplot(gs[0,0], axisbg=bg_color)
    ax1.set_frame_on(True)
    ax1.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax1.yaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax1.spines.values():
	    spine.set_color(fg_color)
    ax1.plot(r1, tmp1_gas ,color=color1, linestyle='-',lw = lw, label= plotting_arguments1['value_label'])
    if data2_dust != []:
        r2 = r_glob2[:,1]*R_ref
        ax1.plot(r2, tmp2_gas, color=color2, linestyle=':', lw = lw, label= plotting_arguments2['value_label'])
    if data3_dust != []:
        r3 = r_glob3[:,1]*R_ref
        ax1.plot(r3, tmp3_gas, color=color3, linestyle='--',lw = lw, label= plotting_arguments3['value_label'])
    if data4_dust != []:
        r4 = r_glob4[:,1]*R_ref
        ax1.plot(r4, tmp4_gas, color=color4, linestyle='-.',lw = lw, label= plotting_arguments4['value_label'])
    if data5_dust != []:
        r5 = r_glob5[:,1]*R_ref
        ax1.plot(r5, tmp5_gas, color=color5, linestyle='-', lw = lw, label= plotting_arguments5['value_label'])
        
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    if data1_gas_norm != []:
        ax1.set_ylabel(r'$\sum_{i=1}^N \Sigma_{g,i}/\Sigma_0 /N $')
    else:
        ax1.set_ylabel(r'$\sum_{i= 1}^N \Sigma_{g,i} /N$')

    ax1.set_ylim([10**(floor(log10(plot_gas_min -0.1* plot_gas_min))), 10**(ceil(log10(plot_gas_max+0.1*plot_gas_max)))])
    ax1.set_xlim([np.floor(r1.min()), np.around(r1.max())])
    ax1.set_xticks([np.floor(r1.min()), R_ref, np.around(r1.max())])    
    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.xaxis.set_ticklabels(['5','','','','','10','','','','','','','','','','','','','','','','','','','','30'],minor=True)
    ax1.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax1.yaxis.set_label_position("left")
    ax1.yaxis.tick_left()

    ax2 = plt.subplot(gs[0,1],axisbg=bg_color)
    ax2.set_frame_on(True)
    for spine in ax2.spines.values():
	    spine.set_color(fg_color)

    ax2.plot(r1, tmp1_dust, color=color1, linestyle='-', lw = lw, label= plotting_arguments1['value_label']) 
    
    if data2_dust != []:   
        ax2.plot(r2, tmp2_dust, color=color2, linestyle=':',  lw = lw, label= plotting_arguments2['value_label'])
    if data3_dust != []:   
        ax2.plot(r3, tmp3_dust, color=color3, linestyle='--', lw = lw, label= plotting_arguments3['value_label'])
    if data4_dust != []:   
        ax2.plot(r4, tmp4_dust, color=color4, linestyle='-.', lw = lw, label= plotting_arguments4['value_label'])
    if data5_dust != []:   
        ax2.plot(r5, tmp5_dust, color=color5, linestyle='-',  lw = lw, label= plotting_arguments5['value_label'])
    
    if data1_gas_norm != []:
        ax2.set_ylabel(r'$\sum_{i=1}^N \Sigma_{d,i}/\Sigma_0 /N $')
    else:
        ax2.set_ylabel(r'$\sum_{i=1}^N \Sigma_{d,i}/N $')
    ax2.set_ylim([10**(floor(log10(plot_dust_min))), 10**(ceil(log10(plot_dust_max)))])

    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlim([np.floor(r1.min()), np.around(r1.max())])
    ax2.set_xticks([np.floor(r1.min()), R_ref, np.around(r1.max())])  
    ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.xaxis.set_ticklabels(['5','','','','','10','','','','','','','','','','','','','','','','','','','','30'],minor=True)
    ax2.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()

    ax1.text(0.98, 0.01, 'Gas surface density',
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax1.transAxes,
        color=fg_color, fontsize=plotting_arguments1['fontsize'])

    ax2.text(0.98, 0.01, 'Dust surface density',
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax2.transAxes,
        color=fg_color, fontsize=plotting_arguments1['fontsize'])

    ax3.text(0.98, 0.01, 'Size',
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax3.transAxes,
        color=fg_color, fontsize=plotting_arguments1['fontsize'])

    ax4.text(0.98, 0.01, 'Stokes number',
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax4.transAxes,
        color=fg_color, fontsize=plotting_arguments1['fontsize'])


    # Put a legend below current axis
    ax4.legend(loc='upper center', bbox_to_anchor=(0, 2.4),
          fancybox=True, shadow=False, ncol=numberforbox)
    plt.savefig('%04d_compare_new_%s.png' %(outerloop, name), facecolor=bg_color, edgecolor=fg_color, dpi=150)
    #plt.savefig('%04d_compare_%s.pdf' %(outerloop, name))
   
    plt.close()





def plot_simple2D(outerloop, plotting_arguments, r, theta, data1, plotrange=[]):
    rcParams.update({'font.size': plotting_arguments['fontsize']}) 
    bg_color    = plotting_arguments['bg_color']
    fg_color    = plotting_arguments['fg_color']
    label_color = plotting_arguments['label_color']
    colormap1   = plotting_arguments['colormap1']
    line_width  = plotting_arguments['line_width']
    name        = plotting_arguments['name']
    nghost 	= plotting_arguments['nghost']
    R_ref 	= plotting_arguments['R_ref']


    plot2D1    = cutforplot(data1, plotting_arguments)


    r_glob                       = cutforplot(r, plotting_arguments)*R_ref
    ytheta                       = np.linspace(0,2*np.pi, len(theta)-2*nghost) 
    radius_matrix , theta_matrix = np.meshgrid(r_glob[:,1],ytheta)

    if plotrange != []:
        vmin1 = plotrange[0]
        vmax1 = plotrange[1] 


    else: 
        vmin1 = plot2D1.min()
        vmax1 = plot2D1.max()

    labelscbar1 = np.linspace(np.round(np.log10(vmin1),1), np.round(np.log10(vmax1),1),num=1+(np.ceil(np.log10(vmax1))- np.floor(np.log10(vmin1))) )


    fig = plt.figure(figsize=(10 , 10),dpi=600, facecolor=bg_color, edgecolor=fg_color)
    gs = GridSpec(2, 1, width_ratios=[1], height_ratios=[15,1])
    gs.update(left=0.01, right=0.99, hspace=0.0, wspace=0.15)

    ax1 = plt.subplot(gs[0,0], projection='polar', axisbg=bg_color)
    ax1.set_frame_on(False)
    ax1.xaxis.set_tick_params(color=fg_color, labelcolor=label_color)
    for spine in ax1.spines.values():
            spine.set_color(fg_color)
    ima = ax1.pcolormesh(theta_matrix, radius_matrix,np.log10(plot2D1.transpose()), 
                    vmin=np.round(np.log10(vmin1),1), 
                    vmax=np.round(np.log10(vmax1),1),cmap=colormap1)
    #ax1.set_rlabel_position(135)
    #rticks = [r_glob.min(), r_glob.max()]
    #rticks = []
    #ax1.set_rgrids(rticks, angle=300)
    ax1.set_yticklabels([])
    ax1.set_theta_zero_location('N')
    thetaticks = []
    ax1.set_thetagrids(thetaticks, frac=1.3)   

    axes1 = plt.subplot(gs[1,0], axisbg=bg_color)
    cba = plt.colorbar(ima,cax=axes1, orientation="horizontal", ticks=labelscbar1)
    cba.ax.yaxis.set_ticks_position('left')
    #cba.set_ticks([-16, -11, -6, -1, 4 ])
    cba.ax.yaxis.set_tick_params(color=fg_color, labelcolor=fg_color, width=0.2)  
    ax1.set_title(r"log$ (c_s)$ ",color=fg_color)    
    #cba.outline.set_visible(False)
    cba.outline.set_linewidth(0.2)
    
    plt.savefig('%04d_%s.png' %(outerloop, name), facecolor=bg_color, edgecolor=fg_color, dpi=600)
    plt.close()
    

