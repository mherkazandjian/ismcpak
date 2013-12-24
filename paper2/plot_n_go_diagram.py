"""
<keywords>
</keywords>
<description>
plot the n-G0 diagram of the paper with the labeled densities, regiems...etc...
</description>
"""
import pylab
import numpy
import os
from mylib.utils import templates

imageSaveDir = '/home/mher/ism/docs/paper02/src/figs'

fw, fh = 3.55, 2.8
text_size = 10.0
clabel_font_size = 8.0

ranges = [[0,6],[0,6]]
cbar_range = [0, 1]
xticks = [0, 1, 2, 3, 4, 5, 6]
yticks = [0, 1, 2, 3, 4, 5, 6]

fig = templates.subplots_grid(1, 1, hspace=0.0, wspace=0.0, sharex=True, sharey=True,
                               fig = {'kwargs':{
                                               'figsize' : (fw, fh),
                                               },
                                      },
                               axs = {
                                      'left':0.2, 'bottom':0.15, 'w':0.7, 'h':0.73
                                      },
                               )

fig.sub_setp('xlim', ranges[0], 'ylim', ranges[1])
fig.set_xticks(xticks, labels=xticks, size=text_size)
fig.set_yticks(yticks, labels=yticks, size=text_size)
fig.set_xlabel('$\log_{10}$ [n$_{gas}$ / cm$^{-3}$]', size=text_size, space=-0.095)
fig.set_ylabel('$\log_{10}$ [G$_0$]'                , size=text_size, space=-0.085)

fig.sub[0,0].plot([0,6], [0,6], 'k--')

fig.sub[0,0].plot([0,3], [3,6], 'k--')

fig.sub[0,0].plot([3,6], [0,3], 'k--')


fig.sub[0,0].text(1, 4.5, r'$G_0 / n = 1000$' , rotation=45, size=8)
fig.sub[0,0].text(3, 3.5, r'$G_0 / n = 1$'    , rotation=45, size=8)
fig.sub[0,0].text(2.8, 1.3, r'$G_0 / n = 0.001$', rotation=40, size=8)


fig.sub[0,0].plot(1.0, 2.0, 'ko' , markersize=3)
fig.sub[0,0].text(1.1, 2.1, 'MA1', size=8)

fig.sub[0,0].plot(2.0, 1.0, 'ko' , markersize=3)
fig.sub[0,0].text(2.1, 1.1, 'MA2', size=8)

fig.sub[0,0].plot(3.0, 3.0, 'ko' , markersize=3)
fig.sub[0,0].text(2.8, 3.3, 'M1' , size=8)

fig.sub[0,0].plot(3.0, 5.0, 'ko' , markersize=3)
fig.sub[0,0].text(3.1, 5.1, 'M2' , size=8)

fig.sub[0,0].plot(5.5, 3.0, 'ko' , markersize=3)
fig.sub[0,0].text(5.5, 3.1, 'M3' , size=8)

fig.sub[0,0].plot(5.5, 5.0, 'ko' , markersize=3)
fig.sub[0,0].text(5.5, 5.1, 'M4' , size=8)


fig.sub[0,0].plot([0,0], [6,6.5], 'k-', clip_on=False)
fig.sub[0,0].text(0.5, 6.3, '  low\ndensity', clip_on=False, size=8)
fig.sub[0,0].plot([2,2], [6,6.5], 'k-', clip_on=False)


fig.sub[0,0].text(3.0, 6.3, 'medium\ndensity', clip_on=False, size=8)
fig.sub[0,0].plot([5,5], [6,6.5], 'k-', clip_on=False)

fig.sub[0,0].text(5.5, 6.3, 'high\ndensity', clip_on=False, size=8)


fig.sub[0,0].text(0.2, 0.2, 'Solar\nNeighbourhood', color='red', size=8)

fig.sub[0,0].text(0.1, 4.5, 'Extermely\nFUV\ndominated', color='red', size=8)


fig.sub[0,0].text(4.0, 5.5, 'Starbursts', color='red', size=8, rotation=-90)

fig.sub[0,0].text(5.1, 5.5, 'Exterme Starbursts', color='red', size=8, rotation=-90)

fig.sub[0,0].text(4.8, 0.5, 'Dark\nClouds', color='red', size=8, rotation=0)

#fig.preview(env='X')
fig.preview(env='aa')


impath = os.path.join(imageSaveDir, 'n_go_diagram.eps')

fig.fig.savefig(impath)

print 'saved image to:\n\t%s' % impath


print 'done'
