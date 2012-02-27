import sys
import pylab as pyl
import numpy as np

t = np.arange(0.0, 1.0, 0.01)
s = np.sin(2*np.pi*t)

fig1 = pyl.figure(1)
ax = fig1.add_subplot(111)
ax.plot(t,s)

def on_move(event):
    # get the x and y pixel coords
    x, y = event.x, event.y

    if event.inaxes:
        ax = event.inaxes  # the axes instance
        print 'data coords', event.xdata, event.ydata

def on_click(event):
    # get the x and y coords, flip y from top to bottom
    x, y = event.x, event.y
    if event.button==1:
        if event.inaxes is not None:
            print 'data coords', event.xdata, event.ydata

binding_id = fig1.canvas.mpl_connect('motion_notify_event', on_move)
binding_id = fig1.canvas.mpl_connect('button_press_event', on_click)


pyl.show()
