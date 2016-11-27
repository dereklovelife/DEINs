# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 21:08:11 2016

@author: d
"""
from math import *
import numpy as np
import numpy.linalg as la
# receive power constrain
R = 0.2
# transmit power constrain
P = 1.0
# BS atenna
N = 3   
# UE atenna
M = 3   
# channel gains
#H = np.array(([0.6,0.65,0.7],[0.63,0.67,0.69],[0.7,0.71,0.63]))
H = np.array(([0.6,0,0],[0,0.58,0],[0,0,0.57]))

# SVD
U,Sigma,Vt = la.svd(H)

Sigma = tuple(Sigma)
h = [i**2 for i in Sigma]



## water filling
mu = h[0]
while 1:
    R2 = R
    ##find lambda
    for i in range(len(h)):
        R2+= log(h[0])-log(h[i])
    lam = exp(R2/len(h))
    
    ##
    if lam/mu-1/h[-1]>0:
        break
    else:
        h.pop()

g = [lam/mu-1/i for i in h]
print g