#! /usr/bin/env python

import math
import numpy as np
import mars_funcs as mf
import matplotlib.pyplot as plt
from pylab import *

CO_flux=np.array((340.75,1813.24,4468.67,5798.58,5612.65,5435.22,4842.61,4625.19,3891.89,3789.22,3141.24,2479.69,2196.16))
yerror=np.array((9.*1.56,49.9,901.88,CO_flux[3]*.3,CO_flux[4]*.3,CO_flux[5]*.3,CO_flux[6]*.3,CO_flux[7]*.3,CO_flux[8]*.3,CO_flux[9]*.3,CO_flux[10]*.3,CO_flux[11]*.3,CO_flux[12]*.3))

file_name='Arp299A_xdr_py'

CO_trans=np.arange(1,CO_flux.size+1)
CO_freq=CO_trans*115.25

z = 0.010300    
CO_freq=CO_freq/(1+z)
CO_flux=CO_flux*3.3284E-23*CO_freq
yerror=yerror*3.3284E-23*CO_freq

#plt.semilogy(CO_flux)
#plt.show()

mods,den,rad=mf.norm_mods(CO_flux)
print mods.shape
rad= np.array(rad)
den= np.array(den)
#plt.semilogy(CO_flux)
#plt.semilogy(mods[-1,])
#plt.semilogy(mods[0,])
#plt.semilogy(mods[150,])
#plt.show()
############################################################
###Start the actual chi_sq reduction part
############################################################
pfac=np.array([0.,0.01,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0,1.1])
radg=np.where(np.array(rad) == 1.6)[0]
last_pdr=radg[0]
chi_sq=[]
p_chi_sq=[]
p_matrix=[]
p_final=[]

for i in range(last_pdr): 
    for j in range(last_pdr):
        for k in range(rad.shape[0]):
            
            for a in range(pfac.shape[0]):
                for b in range(pfac.shape[0]):
                    for c in range(pfac.shape[0]):
                        
                        #print i,j,k,a, b, c
                        #print pfac[a],pfac[b],pfac[c]
                        tot_mod = pfac[c]*mods[k,] + pfac[b]*mods[j,] + pfac[c]*mods[k,]
                        chi_sq1 = np.sum(((CO_flux-tot_mod)**2.)/((13.-9.)*(yerror/CO_flux)**2.))
                        #print chi_sq
                        #plt.semilogy(CO_flux)
                        #plt.semilogy(tot_mod)
                        #plt.semilogy(pfac[c]*mods[k,])
                        #plt.show()
                        p_chi_sq.append([chi_sq1])
                        p_matrix.append([pfac[a],pfac[b],pfac[c]])
            #
                        
            p_chi_sq=np.array(p_chi_sq)
            p_matrix=np.array(p_matrix)
            best=(p_chi_sq.min() == p_chi_sq).nonzero()[0][0]
            #print best
            #print p_matrix.shape
            chi_sq.append([p_chi_sq[best],i,j,k])
            p_final.append(p_matrix[best,0:3])
            #print p_final
            #print chi_sq
            #plt.semilogy(CO_flux)
            #plt.semilogy(p_matrix[best,0]*mods[i,]+p_matrix[best,1]*mods[j,]+p_matrix[best,2]*mods[k,])
            #plt.semilogy(p_matrix[best,2]*mods[k,])
            #plt.show()
            p_matrix=[]
            p_chi_sq=[]
        #end k loop
        
    # end j loop
    
    print i

# end i loop

np.array(chi_sq)
print chi_sq.shape
a=chi_sq[0,:,:,:]
print a.shape
sys.exit()
#best1=(chi_sq.min(0,:,:,:) == p_chi_sq).nonzero()[0][0]

