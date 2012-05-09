import numpy as np
import pylab as pyl
from matplotlib.widgets import Button

freqs = np.arange(2, 20, 3)

ax = pyl.subplot(111)
pyl.subplots_adjust(bottom=0.2)
pyl.axes([-1,7,-1,7])
t = np.arange(0.0, 1.0, 0.001)
s = np.sin(2*np.pi*freqs[0]*t)
l, = pyl.plot(t, s, lw=2)


def next(event):
    pyl.draw()

axnext = pyl.axes([0.81, 0.05, 0.1, 0.075])
bnext = Button(axnext, 'Next')
bnext.on_clicked(next)

pyl.show()