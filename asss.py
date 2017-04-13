# -*- coding: utf-8 -*-
import numpy as np
import numpy.linalg as la
import scipy.stats as st
import matplotlib.pyplot as plt

    
class AntennaSwitch(object):
    def __init__(self, H=None, pmax=10.0, minEnergy=5.0):
        if not H:
            self.numT = 5
            self.numR = 5
            self.H = np.ones((self.numR,self.numT))
            for i in xrange(self.H.shape[0]):
                for j in xrange(self.H.shape[1]):
                    self.H[i, j] = st.rayleigh.rvs()
            self.H = np.mat(self.H)
        else:
            self.H = H
            self.numR, self.numT = self.H.shape
            
        self.cur = [1] * self.numT
        self.cur[0] = 0
        self.pmax = pmax
        self.minEnergy = minEnergy
        # print self.H

    def Iterator(self, shw = True):
        self.result = [0] * (2**self.numT + 1)
        for num in xrange(2**(self.numT) + 1, 2 ** (self.numT+1)):
            s1 = bin(num)[3:]
            self.cur = map(int, list(s1))
            #print self.cur
            self.result[num - 2**(self.numT)] = self.optimize()
        #print self.result
        if shw:
            plt.bar(xrange(2**self.numT + 1), self.result)
            plt.show()
        
        maxvalue = 0
        maxpos = 0
        for i in xrange(len(self.result)):
            if self.result[i] > maxvalue:
                maxvalue = self.result[i]
                maxpos = i
        return (maxvalue, maxpos)

    def optimize(self):
        G = np.mat(np.diag(self.cur))
        curH = self.H * G
        curG = self.H - curH
        
        Hu, sigmaH, Hvt = la.svd(curH)
        Gu, sigmaG, Gvt = la.svd(curG)
        
        sigmaH = sigmaH**2
        sigmaG = sigmaG**2
        
        
        if sigmaG[0] == 0.0 or self.minEnergy/sigmaG[0] > self.pmax:
            return 0.0
        else:
            pleft = self.pmax - self.minEnergy/sigmaG[0]
        PowerAllocation = [0] * self.numT
        PowerAllocation[0] = self.minEnergy/sigmaG[0]
        X1 = np.mat(np.diag(PowerAllocation))
        X1 = Gvt.T * X1 * Gvt        
        #print X1
        #print np.trace(curG*X1*curG.T)
        
        InformationAllocation = [0] * self.numT
        while 1:
            flag = 0
            curSum = 0
            count = 0
            for lamb in sigmaH:
                if lamb == 0:
                    break
                count += 1
                curSum += 1.0/lamb
            Oplamb = count / (pleft+curSum)
            if not Oplamb:
                return 0.0
            ret = 0
            for i in xrange(len(sigmaH)):
                if not sigmaH[i]:
                    flag = 1
                    break 
                if sigmaH[i] < Oplamb:
                    sigmaH[i] = 0
                    break
                ret += np.log(sigmaH[i]/Oplamb)
            if not flag:
                continue
            for i in xrange(len(sigmaH)):
                if not sigmaH[i]:
                    break
                InformationAllocation[i] = 1.0/Oplamb - 1.0/sigmaH[i]
#            X2 = np.mat(np.diag(InformationAllocation))
#            X2 = Hvt.T * X2 *Hvt
#            X = X1 + X2
#            print X
#            print "EH:", np.trace(curG*X*curG.T)
#            print "ID:", np.log(la.det(np.mat(np.eye(self.numT))+curH*X*curH.T)),ret
#            
            return ret
            
    def getH(self):
        return self.H;

    def setMinEnergy(self, n):
        self.minEnergy = n
        
    def RangePlot(self, start = 2, end = 100, n = 20):
        x = np.linspace(start, end, n)
        y = np.ones(n)
        z = np.ones(n)
        for i in xrange(len(x)):
            self.setMinEnergy(x[i])
            y[i], z[i] = self.Iterator(False)
        
        plt.plot(x,y)
        plt.xlabel("Energy Harvested (w)")
        plt.ylabel("Throughput (bps/Hz)")
        plt.grid()
        plt.show()
        plt.plot(x,z,"ro")
        plt.xlabel("Energy Harvested (w)")
        plt.ylabel("Antenna Selection")
        plt.grid()
        plt.show()

if __name__ == "__main__":
    as1 = AntennaSwitch()
    #print as1.getH()
    _,sigma,_ = la.svd(as1.getH())
    print sigma
    as1.RangePlot(2, int(10*sigma[0]**2), 20)
    
