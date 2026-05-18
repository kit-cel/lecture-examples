l# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 14:00:46 2014

@author: jaekel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
   
########################
# main function
########################
def main():
    
    # set time resp. pulse interval and related parameters
    t_min=-5.0
    t_max=5.0
    t_s= 0.01                                 # sample time
    t=np.arange(t_min, t_max+t_s, t_s)
    
    T_rect=2            # width of the rectangular
    t_rect=np.arange(-T_rect/2, T_rect/2+t_s, t_s)
    
    M=len(t)    
    M_rect=len(t_rect)
    
    # rectangular function
    rect = 0*t
    rect[ (M-M_rect)//2 : (M-M_rect)//2+M_rect ] = 1



    # frequency axis and fourier transform
    f_Nyq = 1/(2*t_s)
    delta_f = 1/(t_max-t_min)
    f = np.arange(-f_Nyq, f_Nyq + delta_f, delta_f)    
    
    F_rect = np.fft.fftshift(np.fft.fft(np.roll(rect, (M-1)//2))) 
    
    F_rect = np.fft.fftshift( np.fft.fft( rect ) )
    F_rect = F_rect / np.max(F_rect) * 2.0

    # frequency rect corresponds to time sinc
    Rect_f = 0*f
    
    B = 20 # two-sided bandwidth
    
    Rect_f[ 0 : int(B/delta_f) ] = 1
    Rect_f = np.roll(Rect_f, -int(B/2/delta_f) )
    Rect_f = np.fft.fftshift(Rect_f)


    IF_Rect = np.fft.fftshift( np.fft.ifft( np.fft.fftshift(Rect_f)) )
    IF_Rect = IF_Rect/np.max(IF_Rect)*B
    

    # plotting
    font = {'size'   : 30}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)       
       
    plt.figure(1)
    
    plt.subplot(221)
    plt.plot(t, rect, label='$x(t)$')

    plt.grid(True);    
    plt.xlabel('$t/\mathrm{s}$');    
    #plt.ylabel('$x(t)$')
    plt.legend(loc='upper right')
    plt.axis(xmin=-5, xmax=5, ymin=-0.1, ymax=1.1)
#    plt.annotate('Breite: $T$',
#            xy=(1,0.5),
#            xytext=(2,0.75),
#            arrowprops={'facecolor':'green','shrink':0.05},
#            )



    plt.subplot(222)
    plt.plot(f, F_rect, label='$X(f)$')
    plt.grid(True);    
    plt.xlabel('$f/\mathrm{Hz}$');    
    #plt.ylabel('$X(f)$')
    plt.legend(loc='upper right')    
    #plt.axis(xmin=-25, xmax=25, ymin=-0.5, ymax=2.1)
#    plt.annotate('Nulldurchgang: $1/T$',
#            xy=(1./(T_rect),0.0),
#            xytext=(5,1.0),
#            arrowprops={'facecolor':'green','shrink':0.05},
#            )
            
    plt.subplot(223)
    plt.plot(t, IF_Rect, label='$x(t)$')
    plt.grid(True);    plt.xlabel('$t/\mathrm{s}$');    plt.ylabel('$x(t)$')
    plt.axis(xmin=-2, xmax=2, ymin=-0.25*20, ymax=1.1*20) 
    plt.legend(loc='upper right')     

    plt.subplot(224)
    plt.plot(f, Rect_f, label='$X(f)$')
    plt.grid(True);    plt.xlabel('$f/\mathrm{Hz}$');    plt.ylabel('$X(f)$')
    plt.axis(xmin=-25, xmax=25, ymin=-0.1, ymax=1.1)
    plt.legend(loc='upper right')        

   
  
    
    plt.show()
    

  
########################
# make it executable
########################
if __name__ == "__main__":
    main()