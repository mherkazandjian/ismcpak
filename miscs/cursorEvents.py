asdasdasd
# MAKE SURE THE MESH DATA CORRESPONDS TO THE ORIGINAL ONES!!!!!!!!
# USE ORIGINAL INDICIES

def onMove(event):
        # get the x and y coords, flip y from top to bottom
    b      = event.button 
    x, y   = event.x, event.y
    xd, yd = event.xdata, event.ydata
        
    if event.button==1:
        if event.inaxes is not None:
            print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % (b, x, y, xd, yd)

            
            l2Distance  = np.sqrt( (xd - lG0 )**2 + (yd - lnGas)**2 )
            rMin = min(l2Distance)
            indMin = l2Distance.argmin()

            """
            mDummy.clearLinesInSubPlots()
            mDummy.data = arxv.meshes[indMin][0]
            mDummy.setFigure(fig, axs)
            mDummy.plot(eSpcs)
            fig.show()
            """

            plt.figure(1)
            plt.hold(False)
            plt.plot( lG0, lnGas, 'r+' )
            plt.hold(True)
            plt.plot( lG0[indMin], lnGas[indMin], 'g+', markersize=20)
            plt.show()            

 
def onClick(event):
    # get the x and y coords, flip y from top to bottom
    b      = event.button 
    x, y   = event.x, event.y
    xd, yd = event.xdata, event.ydata
        
    if event.button==1:
        if event.inaxes is not None:
            
            l2Distance  = np.sqrt( (xd - lG0 )**2 + (yd - lnGas)**2 )
            rMin = min(l2Distance)
            indMin = l2Distance.argmin()

            #plt.figure(2)
            mDummy.clearLinesInSubPlots()
            mDummy.data = arxv.meshes[indMin][0]
            #mDummy.setFigure(fig, axs)
            #mDummy.plot(eSpcs)
            #plt.show()
            
            #plt.figure(1)
            #plt.hold(False)
            #plt.plot( lG0, lnGas, 'r+' )
            #plt.hold(True)
            #plt.plot( lG0[indMin], lnGas[indMin], 'g+', markersize=20)
            #plt.show()

                        
cid = fig1.canvas.mpl_connect('button_press_event' , onClick)
cid = fig1.canvas.mpl_connect('motion_notify_event', onMove)
fig1.show()
