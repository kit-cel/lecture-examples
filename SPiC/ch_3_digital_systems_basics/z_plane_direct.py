# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 08:34:10 2014

@author: jaekel
"""

########################
######## importing
########################
  
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


########################
# main function
########################
def main():
    
    # define parameters
    zeros = np.array([0, 0, 0])
    poles = np.array([.5, .25+.8j, .25-.8j])

    assert len(zeros)==len(poles)
    
    Omega = np.linspace(-np.pi, np.pi, 512)    
    H = np.zeros( len(Omega), dtype=complex )    
    
    for O in Omega:
        H[ np.where(Omega==O) ] = np.prod( np.exp(1j*O)*np.ones(len(zeros)) - zeros) / np.prod( np.exp(1j*O)*np.ones(len(poles)) - poles)
    
    
    # plotting
    font = {'size'   : 26}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)  

    
    plt.figure(2)    
    plt.subplot(121)
    plt.plot(Omega, abs(H))
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    plt.xlabel('$\Omega$')
    plt.ylabel('$|H(\Omega)|$')

    plt.subplot(122)    
    plt.plot(Omega, np.angle(H))
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    plt.xlabel('$\Omega$')
    plt.ylabel('$\phi(\Omega)$')    

    plt.show()
    

    # that's it folks!
    print('Done!')


########################
# make it executable
########################
if __name__ == "__main__":
    main()
