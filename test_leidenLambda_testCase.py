# test with a 2 and 3 level system


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
#### NEED TO MODIFY THIS WITH THE NEW COLLISION FORMAT OF moldata.py
####-----------------------------------------------------------------

#parameters
Tkin  = 500.0    # kinetic temperature of the gas (K)
Tcmb  = 2.73     # temperature of background radiation (K)
nc    = 1e12     # collider number density cm^-3
# some constants
#------------------------------
hPlank = 6.63e-27       # erg.s 
cLight = 29979245800.0  # cm.s^-1
kBoltz = 1.38e-16       # erg.K^-1 
ev2erg = 1.602e-12      #erg
#------------------------------

#reading the whole database of line info of species from LAMBDA
lambdaPath = '/home/mher/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles'
reader = molData.reader(dirPath = lambdaPath, species = 'HCN')

spec = reader.get_specie(specStr='HCN', inPath='hcn.dat')


#number of energy levels
n     = spec.nlevels
    
# converting the energy units to K 
for idx, level in enumerate(spec.levels):
    level['E'] *= hPlank*cLight/kBoltz # energies in K
levels = spec.levels
transRad = spec.transRad
transColl = spec.transColl['H2']['trans']

print '----------------------------levels----------------------------'
print 'idx         E      g'
print '--------------------------------------------------------------'
for l in levels:
    print '%03d     %.2e   %d ' % (l['n'], l['E'], l['g'])  

print '------------------radiative transitions--------------------------'
print 'idx   u    l    A          nu        E'
print '-----------------------------------------------------------------'
for t in transRad:
    print '%d     %d    %d  %.2e  %.2e %.2e' % (t['n'], t['u'], t['l'], t['A'], t['nu'], t['E'])
      
print '------------------collisional transitions-----------------------'
print 'idx       u      l    rc(T=10) rc(T=15) rc(T=20)'
print '-----------------------------------------------------------------'
for t in transColl:
    print '%05d     %03d    %03d  %.2e %.2e %.2e' % (t['n'], t['u'], t['l'], t['rc'](10.0), t['rc'](15.0), t['rc'](20.0))  

###########################################################
# constructing the matrix 
###########################################################

#a utility function which prints a matrix
def print_matrix( m ):
    print '     ',
    for i, row in enumerate(m):
        print '  %02d   ' % i,
    print 
    for i in np.arange(100):
        print '-',
    print 
    
    for i, row_i in enumerate(m):
        print '%02d   ' % i,
        for mij in row_i:
            print '%.1e' % mij,
        print
     
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
    
    #print_matrix(K)        
    
    return K


def fill_AP_matrix():
    """fill the Aij matrix from the lambda radiative transitions database"""
    AP = np.zeros( (n, n), dtype = np.float64)
    
    for trans in transRad:
        u  = trans['u']; l = trans['l']
        dE = np.abs( levels[u]['E'] - levels[l]['E'] ) # energy in K
        
        nu = dE*kBoltz / hPlank # freq in Hz
        
        AP[u, l] = (1.0 + ng(hPlank, nu, kBoltz, Tcmb))*trans['A']

    #print_matrix(AP)        
    
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

    #print_matrix(ABS)        
    
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
dndt = np.zeros((n,1))

# solving directly
#replacing the first row with the conservation equation
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

for i, level in enumerate(levels):
    print 'i = %02d, level = %02d, pop_dens = %e' % (i, level['n'], x[i]) 
    
print 'done'