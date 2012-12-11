#this is a test file for leidenLambda.py.
import sys, os
if 'particle3' in os.uname():
    import matplotlib
    matplotlib.use('Qt4Agg')
import pylab as pyl

from leidenLambda import molData 
import numpy as np
import re
from ismUtils import planckOccupation as ng

restore = False
#reading the whole database of line info of species from LAMBDA
lambdaPath = '/home/mher/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles'
reader = molData.reader(dirPath = lambdaPath)

Tkin_min    = 15.0
Tkin_max    = 300.0
TkinSamples = 40.0 
nH2         = [1e2, 1e4, 1e6]

Tcmb     = 2.73  # temperature of background radiation (K)
# some constants
#------------------------------
hPlank = 6.63e-27    # erg.s 
cLight = 29979245800.0 # cm.s^-1
kBoltz = 1.38e-16    # erg.K^-1 
ev2erg = 1.602e-12   #erg
#------------------------------

#selecting the one holding the info for p-NH3
for specDict in reader.speciesInfo:
    if 'p-NH3' in specDict['specStr']:
        pNH3 = specDict
        print 'found p-NH3 info, filePath : %s' % pNH3['path']
        break
    
# converting the energy units to K 
for idx, level in enumerate(pNH3['levels']):
    level['E'] *= hPlank*cLight/kBoltz # energies in K

#############################################################################
# compute the matrix of transition rates
def computeRateMatrix(pNH3, Tkin, nc):
    n = pNH3['nlevels']
         
    levels = pNH3['levels']
    transRad = pNH3['transRad']
    transColl = pNH3['transColl']['partner0']['trans']
      
    ###########################################################
    # constructing the matrix 
    ###########################################################
    def fill_K_matrix():
        """fill the kij matrix from the lambda collisional transition database"""
        #n  = 5
        K = np.zeros( (n, n), dtype = np.float64)
        for trans in transColl:
            u  = trans['u']; l = trans['l']
            gu = levels[u]['g']; gl = levels[l]['g']
            
            dE = np.abs(levels[u]['E'] - levels[l]['E']) # difference in the enrergy in K
            
            K[u, l] = trans['rc'](Tkin)
            K[l, u] = (np.float64(gu) / np.float64(gl)) * K[u,l] * np.exp(-dE / Tkin)
        
        return K
    
    
    def fill_AP_matrix():
        """fill the (A prime)_ij matrix from the lambda radiative transitions database"""
        AP = np.zeros( (n, n), dtype = np.float64)
        
        for trans in transRad:
            u  = trans['u']; l = trans['l']
            dE = np.abs( levels[u]['E'] - levels[l]['E'] ) # energy in K
            
            nu = dE*kBoltz / hPlank # freq in Hz
            
            AP[u, l] = (1.0 + ng(hPlank, nu, kBoltz, Tcmb))*trans['A']
    
        return AP
    
    def fill_ABS_matrix():
        """fill the Aij matrix for absorbtion transitions from the lambda radiative transitions database"""
        ABS = np.zeros( (n, n), dtype = np.float64)
        
        for trans in transRad:
            u  = trans['u']; l = trans['l']
            dE = np.abs( levels[u]['E'] - levels[l]['E'] ) # energy in K
            gu = levels[u]['g']; gl = levels[l]['g'] 
            
            nu = dE*kBoltz / hPlank # freq in Hz
            
            ABS[u, l] = (np.float64(gu)/np.float64(gl))*ng(hPlank, nu, kBoltz, Tcmb)*trans['A']
    
        return ABS
    
    def fill_E_matrix():
        """This is a matrix full of ones (it is a utility matrix"""
        return np.ones((n,n))
    
    K = fill_K_matrix()
    AP = fill_AP_matrix()
    ABS = fill_ABS_matrix()
    E = fill_E_matrix()
    
    F = nc * K + AP + ABS.T
    diag = np.eye(n)*np.dot(F,E)[:,1]
    offdiag = F.T
    
    full = -diag + offdiag

    return full
#############################################################################
# solve for the equlibrium population densities
def solveEquilibrium(pNH3, full):
    n = pNH3['nlevels']

    # solving directly
    #replacing the first row with the conservation equation
    dndt = np.zeros((n,1))
    full[0,:] = 1.0
    dndt[0]   = 1.0
    
    A = full
    b = dndt
    #solving the system A.x = b
    #--------------------------
    #before solving, we will devide each row by the diagonal
    for i in np.arange(n):
        A[i,:] = A[i,:]/A[i,i]
    x = np.linalg.solve(A, b)
    
    #print x.T
    
    # the fractional population density
    f = x
    return f    
#############################################################################
def computeLuminosity(pNH3):
    n = pNH3['nlevels']
    levels = pNH3['levels']
    transRad = pNH3['transRad']

    ## computing the total luminosity##
    #----------------------------------
    def fill_R_matrix():
        """fill the L_ij matrix from the lambda radiative transitions database"""
        R = np.zeros( (n, n), dtype = np.float64)
        
        for trans in transRad:
            u  = trans['u']; l = trans['l']
            dE = np.abs( levels[u]['E'] - levels[l]['E'] ) # energy in K
            
            nu = dE*kBoltz / hPlank # freq in Hz
            
            R[u, l] = trans['A']*hPlank*nu
    
        return R
    #------------------------------------
    
    R = fill_R_matrix()
    
    return R
#############################################################################

####### plotting the line intensities per molecule as a function of Tkin#####
for nc in nH2:
    
    Tkins = np.linspace(Tkin_min, Tkin_max, TkinSamples)
    Llines = {'11':[], '22':[], '44':[]}
    
    for Tkin in Tkins:
        
        full = computeRateMatrix(pNH3, Tkin, nc)
        f    = solveEquilibrium(pNH3, full)
        R    = computeLuminosity(pNH3)
        
        u11 = 1; l11 = 0;    
        L11 = R[u11][l11]*f[u11]
        
        u22 = 3; l22 = 2;    
        L22 = R[u22][l22]*f[u22]
        
        u44 = 11; l44 = 10;    
        L44 = R[u44][l44]*f[u44]
    
        Llines['11'].append(L11)
        Llines['22'].append(L22)
        Llines['44'].append(L44)
        
    if nc == 100.0:
        colorStr = 'g'
    if nc == 10000.0:
        colorStr = 'r'
    if nc == 1000000.0:
        colorStr = 'b'
        
    pyl.plot(Tkins, Llines['11'], colorStr)
    pyl.hold(True)
    pyl.plot(Tkins, Llines['22'], colorStr+'--')
    pyl.plot(Tkins, Llines['44'], colorStr+'-.')
    
pyl.xscale('log')    
pyl.yscale('log')
pyl.axis([15,300,1e-25,1e-22])
pyl.show()

"""
################ plotting the line ratios as a function of Tkin###############
for nc in nH2:
    
    Tkins = np.linspace(Tkin_min, Tkin_max, TkinSamples)
    Llines = {'11':[], '22':[], '44':[]}
    
    for Tkin in Tkins:
        
        full = computeRateMatrix(pNH3, Tkin, nc)
        f    = solveEquilibrium(pNH3, full)
        R    = computeLuminosity(pNH3)
        
        u11 = 1; l11 = 0;    
        L11 = R[u11][l11]*f[u11]
        
        u22 = 3; l22 = 2;    
        L22 = R[u22][l22]*f[u22]
        
        u44 = 11; l44 = 10;    
        L44 = R[u44][l44]*f[u44]
    
        Llines['11'].append(L11)
        Llines['22'].append(L22)
        Llines['44'].append(L44)
        
    if nc == 100.0:
        colorStr = 'g'
    if nc == 10000.0:
        colorStr = 'r'
    if nc == 1000000.0:
        colorStr = 'b'
        
    pyl.plot(Tkins, np.array(Llines['22'])/np.array(Llines['11']), colorStr)
    pyl.hold(True)
    pyl.plot(Tkins, np.array(Llines['44'])/np.array(Llines['22']), colorStr+'--')
    pyl.plot(Tkins, np.array(Llines['44'])/np.array(Llines['11']), colorStr+'-.')
    
pyl.xscale('log')    
pyl.yscale('log')
pyl.axis([15, 300, 0.01, 1.0])
pyl.show()    
"""

print 'done'