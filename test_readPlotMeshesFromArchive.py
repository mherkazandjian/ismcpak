from mesh import *
from meshUtils import *
from enumSpecies import *
import numpy as np
from time import *
import pylab as pyl
 
nMove = 0
indMin = 0

runDirPath    =  '/home/mher/ism/runs/oneSided/uniformSweep2-z-2/'
abunSpecFname = 'data/species.inp'
lgammaMechSec = -22.0

# reading the archive 
t0 = time()
arxv = meshArxv(  )
arxv.readDb( runDirPath )
arxv.checkIntegrity()
print 'time reading %f' % (time() - t0)

specFileUsed = enumSpecies(abunSpecFname, offsetFromZero=1 )
eSpcs = specFileUsed.genPythonClass(offset=0)

lG0All    = np.log10(arxv.infoAll['parms'][:,0])
lnGasAll  = np.log10(arxv.infoAll['parms'][:,1])
lgmAll    = np.log10(arxv.infoAll['parms'][:,2])
mDummy = mesh()
nMove  = 0

# check if the information about the meshes in the header of the 
# archive and the same information in the database are the same

indsThisSec = ( np.fabs(lgmAll - lgammaMechSec) < 1e-6 ).nonzero()


fig1 = pyl.figure(1, figsize=(6, 6))
axs1 = fig1.add_subplot(111)
axs1.plot( lG0All[indsThisSec], lnGasAll[indsThisSec], 'bo' )

fig2, axs2 = pyl.subplots(2, 2, sharex=True, sharey=False, figsize=(8,8))


def onB1DownMove(event):
    global nMove
    nMove += 1
    
    # get the x and y coords, flip y from top to bottom
    b      = event.button 
    x, y   = event.x, event.y
    xd, yd = event.xdata, event.ydata
        
    if event.button==1:
        if event.inaxes is not None:
            
            l2Distance  = np.sqrt( (xd - lG0 )**2 + (yd - lnGas)**2 )
            rMin = min(l2Distance)
            indMin = l2Distance.argmin()

            print xd, yd, nMove

            pyl.figure(1)
            pyl.subplot(111)
            pyl.hold(False)
            pyl.plot( lG0, lnGas, 'bo' )
            pyl.hold(True)
            pyl.plot( lG0[indMin], lnGas[indMin], 'ro')
            pyl.title('$\log_{10} G\_0 = $ %4.2f $\log_{10} n_{gas} = $ %4.2f  $\log_{10} \Gamma_{mech} = $  %5.2f\n ' % (xd, yd, lgammaMechSec) )
            pyl.draw()

            pyl.figure(2)
            pyl.clf()
            
            fig2, axs2 = pyl.subplots(2, 2, sharex=True, sharey=False)
            #pyl.subplots(2, 2, sharex=True, sharey=False)
            
            axs2[0,0]   .plot([1,2,3,4])
            """
            fig2.add_subplot(221)
            fig2.add_subplot(222)
            fig2.add_subplot(223)
            fig2.add_subplot(224)
            """
            
            """
            pyl.title('asdasdasd')
            mDummy.setData( arxv.meshes[indMin][0] )
            mDummy.clearLinesInSubPlots()
            mDummy.setFigure(fig2, axs2)
            mDummy.plot(eSpcs)
            """
            pyl.draw()

def onB1Down(event):
    
    # get the x and y coords, flip y from top to bottom
    b      = event.button 
    x, y   = event.x, event.y
    xd, yd = event.xdata, event.ydata
        
    if event.button==1:
        if event.inaxes is not None:
            
            l2Distance  = np.sqrt( (xd - lG0All )**2 + (yd - lnGasAll)**2 + (lgammaMechSec - lgmAll)**2 )
            rMin = min(l2Distance)
            indMin = l2Distance.argmin()

            #print xd, yd, nMove

            pyl.figure(1)
            pyl.subplot(111)
            pyl.hold(False)
            pyl.plot( lG0All[indsThisSec], lnGasAll[indsThisSec], 'bo' )
            pyl.hold(True)
            pyl.plot( lG0All[indMin], lnGasAll[indMin], 'ro')
            pyl.title('$\log_{10} G_0 = $ %4.2f $\log_{10} n_{gas} = $ %4.2f  $\log_{10} \Gamma_{mech} = $  %5.2f\n ' % (xd, yd, lgammaMechSec) )
            pyl.xlabel('$log_{10} G_0$')
            pyl.ylabel('$log_{10} n_{Gas}$')
            pyl.draw()

            pyl.figure(2)            
            mDummy.setData( arxv.meshes[indMin][0] )
            mDummy.setFigure(fig2, axs2)
            mDummy.plot(eSpcs)

            pyl.draw()

#cid = fig1.canvas.mpl_connect('motion_notify_event', onB1DownMove)
cid = fig1.canvas.mpl_connect('button_press_event', onB1Down)

pyl.figure(1)
pyl.show()


print 'done'