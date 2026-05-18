# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 10:10:34 2013

@author: jaekel
"""

########################
# importing
########################
import numpy as np
import scipy as sp
from numpy import random

import pylab as pl

# define dimensions of signal and cs
n=1000
m=10


# misc. parameters
max_degree=5
degree=30#random.randint(1, max_degree)

t_max=1
t_step=t_max/float(n)
t=np.arange(0, t_max, t_step)

def main():
    
    # construct random signal by randomly choosing a polynmial degree and coefficients 
    coeff=5.0*random.randn(degree)+.5
    
    f=np.zeros(n)
    k=0
    
    while k<degree:
        f=f+coeff[k]*(t**k)
        k+=1
        
    pl.plot(t, f)
    pl.show()
    



########################
# make it executable
########################
if __name__ == "__main__":
    main()