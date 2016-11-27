# -*- coding: utf-8 -*-
"""
Created on Fri Nov 04 11:44:03 2016
This script is to numberical simulate papaer 
"MIMO Broadcasting for Simultaneous
Wireless Information and Power Transfer"

@author: d
"""
import random
from math import *
import numpy as np
import numpy.linalg as la

## BS antennas
n = 2       

## UE antennas
m = 2

## Total power restrains
P_sum = 1

def Channal(n,m):
    """
    Establish Channal Model 
    Input int n, int m
    Return m*n Channal gains matrix H
    """
    
    H = np.ones((n,m))
    H[(0,0)] = 0.8
    H[(0,1)] = H[(1,0)] = 0
    H[(1,1)] = 0.75
    print H
    return H
    
def maxQ():
    H = Channal(n,m)
    U,sigma,Vt = la.svd(H)
    print "Optimal MAX Q:",sigma[0]**2*P_sum
    S1 = np.matrix(np.zeros((2,2)))
    
    S1[(0,0)] = P_sum
    print S1
    S = Vt.transpose()*S1*Vt
    print "Optimal S:",S
    print "Q:",np.trace(H*S*H.transpose())
    print "R:",log(la.det(np.identity(m)+H*S*H.transpose()))
    
def maxR():
    H = Channal(n,m)
    U,sigma,Vt = la.svd(H)
    print "sigma:",sigma
    ## optimal list 
    table = sigma
    S1 = np.matrix(np.zeros((n,n)))
    ret = [-1]*n
    while 1:
        flag =1
        count = 0 
        sum1 = 0
        
        for i in range(n):
            if ret[i] != 0:
                if table[i]!=0:
                    sum1+=1.0/table[i]**2
                    count+=1
                else:
                    ret[i] = 0
                    
        k = (P_sum+sum1)/count
        
        for i in range(n):
            if ret[i] != 0:
                ret[i] = k-1.0/table[i]**2
                if ret[i] <=0:
                    ret[i] = 0
                    flag = 0
        if flag == 1:
            break
    for i in range(n):
        S1[(i,i)] = ret[i]
    ##print "S1:",S1
    S = Vt.transpose()*S1*Vt
    print "S:",S
    print "optimal R:",log(la.det(np.identity(m)+H*S*H.transpose()))
    print "Q:",np.trace(H*S*H.transpose())
    
maxQ()
maxR()
    
    