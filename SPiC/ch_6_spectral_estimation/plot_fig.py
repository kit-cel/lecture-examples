# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 09:02:04 2014

@author: jaekel
"""

########################
# plotting data given in txt format
########################




# importing relevant stuff
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def main():

    # define file name including path here
    file_dir = './'
    file_name = 'results'   
    
    # read from file
#    file_obj = open( file_dir + file_name, 'r')
#    
#    while 1:    
#        line = file_obj.readline()
#        if len(line)==0:
#            break
#        i = 0
#        if line[0]!='-' and i<2:
#            print line.split(';;;')
#            i += 1
    
    data = np.load(file_name, x)#N_acf, y=acf_unbiased)
    print x
    
    ######################################      
    # plotting
    ######################################  
    if 0:
        font = {'size': 26}
        matplotlib.rc('font', **font)
        matplotlib.rc('text', usetex=True)      
        
        # plot acfs
        plt.figure(1)    
        plt.subplot(111)
        plt.plot(x_1, y_1)      

      
        
        plt.show()
        



    
########################
# make it executable
########################
if __name__ == "__main__":
    main()