# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 08:36:12 2014

@author: jaekel
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
   
########################
# main function
########################
def main():
    
    # set time resp. pulse interval and related parameters
    t_min=-0.0
    t_max=10.0
    t_s= 0.01                                 # sample time
    t=np.arange(t_min, t_max+t_s, t_s)
    
    f0=1 
    
    t_sample_1=.1    
    t_sampled_1=np.arange(t_min, t_max+t_sample_1, t_sample_1)
    
    t_sample_2=.4
    t_sampled_2=np.arange(t_min, t_max+t_sample_2, t_sample_2)
    
    # frequency regime
    #f=np.arange( -1/(2*t_s), 1/(2*t_s)+1/(t_max-t_min), 1/(t_max-t_min))
    
    # original signal    
    x=np.sin(2*np.pi*f0*t)
    x_sampled_1=np.sin(2*np.pi*f0*t_sampled_1)
    x_sampled_2=np.sin(2*np.pi*f0*t_sampled_2)    

  

    # plotting
    font = {'size'   : 36}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)
       

    plt.figure(1)
    
    plt.subplot(121)
    plt.plot(t, x, label='$x(t)$')
    plt.plot(t_sampled_1, x_sampled_1,'ro', markersize=10)
    plt.grid(True); plt.xlabel('$t/\mathrm{s}$');  #plt.ylabel('$x(t)$')
    plt.axis(xmin=0, xmax=1.5, ymin=-1.1, ymax=1.1)
    plt.legend(loc='upper right')
    plt.annotate('T=0.1',
            xy=(.6,0.75),
            xytext=(.4,.75),
            )     
    

    plt.subplot(122)
    
    x_recon_1=0*t
    
    for k in np.arange(0, len(t_sampled_1)):
        argument = np.pi/t_sample_1*(t-k*t_sample_1)

        recon_sum_part_1 = np.sin( argument ) / argument
        positions_nan_1 = np.isnan(recon_sum_part_1)
        recon_sum_part_1[ positions_nan_1 ] = 1
        
        x_recon_1 += x_sampled_1[k] * recon_sum_part_1
        
        plt.plot(t, x_sampled_1[k] * recon_sum_part_1,'--')
    
    plt.plot(t,x_recon_1, label='$x(t)$')
    plt.grid(True); plt.xlabel('$t/\mathrm{s}$'); # plt.ylabel('$x(t)$')
    plt.axis(xmin=0, xmax=.5, ymin=-.5, ymax=1.1)
    plt.legend(loc='upper right')    
    plt.annotate('Rekonstruiertes Signal',
            xy=(.14,0.8),
            xytext=(.2,-.3),
            arrowprops={'facecolor':'green','shrink':0.05},
            )
    plt.annotate('sinc-Kerne',
            xy=(.06,-.2),
            xytext=(.1,-.4),
            arrowprops={'facecolor':'green','shrink':0.05},
            )            
    
    

    plt.figure(2)
    
    plt.subplot(121)
    plt.plot(t, x, label='$x(t)$')
    plt.plot( t_sampled_2, x_sampled_2,'ro', markersize=10)
    plt.grid(True); 
    plt.xlabel('$t/\mathrm{s}$');  
    #plt.ylabel('$x(t)$')
    plt.axis(xmin=0, xmax=1.5, ymin=-1.1, ymax=1.1)
    plt.annotate('T=0.4',
            xy=(.5,0.75),
            xytext=(.4,.75),
            )     
    plt.legend(loc='upper right')
    

    plt.subplot(122)
    
    x_recon_2=0*t
    for k in np.arange(0, len(t_sampled_2)):
        argument = np.pi/t_sample_2*(t-k*t_sample_2)

        recon_sum_part_2 = np.sin( argument ) / argument
        positions_nan_2 = np.isnan(recon_sum_part_2)
        recon_sum_part_2[ positions_nan_2 ] = 1
        
        x_recon_2 += x_sampled_2[k] * recon_sum_part_2
        
        plt.plot(t, x_sampled_2[k] * recon_sum_part_2,'--')
    
    plt.plot(t,x_recon_2, label='$x(t)$')
    plt.grid(True); 
    plt.xlabel('$t/\mathrm{s}$');  
    #plt.ylabel('$x(t)$')
    plt.axis(xmin=0, xmax=.5, ymin=-.5, ymax=1.1)
    plt.annotate('Rekonstruiertes Signal',
            xy=(.1,0.35),
            xytext=(.2,-.3),
            arrowprops={'facecolor':'green','shrink':0.05},
            )
    plt.annotate('sinc-Kerne',
            xy=(.1,-.05),
            xytext=(.15,-.4),
            arrowprops={'facecolor':'green','shrink':0.05},
            )   
    plt.legend(loc='upper right')
    
    
    plt.show()
    

  
########################
# make it executable
########################
if __name__ == "__main__":
    main()