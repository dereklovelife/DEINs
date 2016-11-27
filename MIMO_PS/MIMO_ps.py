# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 19:18:35 2016

@author: d

This simulate point-to-point MIMO DEIN with PS approach.



"""
from math import *
import numpy as np
import numpy.linalg as la

# receive power constrain
Q = 0.2
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

P1 = P-Q/h[0]

## water filling
while 1:
    P2 = P1
    for i in range(len(h)):
        P2+=1/h[i]
    mu = len(h)/P2
    if 1/mu-1/h[-1]>0:
        break
    else:
        h.pop()

g = [1/mu-1/i for i in h]
print g




