# -*- coding: utf-8 -*-
"""
Created on Fri Nov 04 15:29:37 2016

@author: d
"""

from math import *
import numpy as np
import numpy.linalg as la

## BS antennas
n = 2       

## UE antennas
m = 2

## Total power restrains
P_sum = 1

Q = 0.60

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
    #print H
    return H

def findu(h,n,P,Q):
    u=1.0
    du = 0.005
    while 1:
        p = findp(u,h,n,Q)
        sbg = P-sum(p)
        if abs(sbg)<0.005:
            break
        u-=du*sbg
    return p
    
def findp(u,h,n,Q):

    sum2 = 0
    p =[0]*n
    hmax = max(h)
    lamb_max = u/hmax
    lamb_min = 0
    while (lamb_max-lamb_min>0.0005):
        lamb = (lamb_max+lamb_min)*0.5
        sum1 = 0
        for i in range(n):
            sum1 += h[i]/(u-lamb*h[i])/log(2)
        if sum1<Q+n:
            lamb_min = lamb
        else:
            lamb_max =lamb
    
    for i in range(n):
        p[i] = 1/(u-lamb*h[i])/log(2) - 1.0/h[i]
        if p[i]<0:
            print "gg"
    return p
    
def main():
    H = Channal(n,m)
    U,sigma,Vt = la.svd(H)
    h = [sigma[i]**2 for i in range(n)]
    p = findu(h,n,P_sum,Q)
    S1 = np.matrix(np.zeros((n,m)))
    for i in range(n):
        S1[(i,i)] = p[i]
    S = Vt.transpose()*S1*Vt
    
    print "Q:",np.trace(H*S*H.transpose())
    print "Q~:",np.trace(H*U*S*U.transpose()*H.transpose())
    print "R:",log(la.det(np.identity(n)+H*S*H.transpose()))
    print "c(y):",H*S*H.transpose()
    print "S:",S
    
main()