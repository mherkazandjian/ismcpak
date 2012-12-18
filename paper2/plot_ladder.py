# plot the ladder of transitions (tested for CO data) emission and optical depth
#------------------------------------------------------------------------------------
import numpy as np
import pickle, os
if 'particle3' in os.uname():
    import matplotlib
    matplotlib.use('Qt4Agg')
import pylab as pyl

dataDir = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/CO/ladder/n-4-g0-4/'
files  = [ '1e-3', '1e-2', '5e-2', '1e-1'] #gmech fractions of surface heating
imageSavePath = '/home/mher/ism/docs/paper02/src/figs/CO-ladder.eps'
#====================================================================================== 

fig, axs = pyl.subplots(1, 2, figsize = (12, 6))

pyl.subplot(121)
data = np.loadtxt(dataDir + 'gm-' + files[0] + '/intensities.out'); 
x = data[:,0]; y = data[:,1]
p1, = pyl.semilogy(x, y, 'r')
#-----------------------
data = np.loadtxt(dataDir + 'gm-' + files[1] + '/intensities.out')
x = data[:,0]; y = data[:,1]
p2, = pyl.semilogy(x, y, 'b')
#-----------------------
data = np.loadtxt(dataDir + 'gm-' + files[2] + '/intensities.out')
x = data[:,0]; y = data[:,1]
p3, = pyl.semilogy(x, y, 'g')
#-----------------------
data = np.loadtxt(dataDir + 'gm-' + files[3] + '/intensities.out')
x = data[:,0]; y = data[:,1]
p4, = pyl.semilogy(x, y, 'c')

pyl.axis([0, 20, 1e-10, 1])
pyl.legend([p1, p2, p3, p4], ['0%','0.1%','1%','10%'])
pyl.xlabel('J')
pyl.ylabel('$\log_{10} [ flux / erg.cm^{-2}.s^{-1} ]$')
#####################################################################
#####################################################################
#####################################################################
pyl.subplot(122)

data = np.loadtxt(dataDir + 'gm-' + files[0] + '/tau.out'); 
x = data[:,0]; y = data[:,1]
p1, = pyl.plot(x, y, 'r')
#-----------------------
data = np.loadtxt(dataDir + 'gm-' + files[1] + '/tau.out')
x = data[:,0]; y = data[:,1]
p2, = pyl.plot(x, y, 'b')
#-----------------------
data = np.loadtxt(dataDir + 'gm-' + files[2] + '/tau.out')
x = data[:,0]; y = data[:,1]
p3, = pyl.plot(x, y, 'g')
#-----------------------
data = np.loadtxt(dataDir + 'gm-' + files[3] + '/tau.out')
x = data[:,0]; y = data[:,1]
p4, = pyl.plot(x, y, 'c')

pyl.axis([0, 20, 0, 500])
pyl.legend([p1, p2, p3, p4], ['0%','0.1%','1%','10%'])
pyl.xlabel('J')
pyl.ylabel('$tau$')




fig.savefig(imageSavePath)
print 'wrote image : %s' % imageSavePath
pyl.show()

print 'done'
