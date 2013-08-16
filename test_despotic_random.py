from despotic import cloud

from numpy import random

import time


mycloud = cloud()

t0 = time.time()

nTrials = 0

while nTrials < 100: 
    mycloud.nH = 1.0e4*(1.0 + random.rand()*0.1)
    mycloud.colDen = 1.0e22*(1.0 + random.rand()*0.1)
    mycloud.sigmaNT = 4.0e5*(1.0 + random.rand()*0.1)
    mycloud.Tg = 300.0*(1.0 + random.rand()*0.1)
    mycloud.comp.xoH2 = 0.1*(1.0 + random.rand()*0.1)
    mycloud.comp.xpH2 = 0.4*(1.0 + random.rand()*0.1)
    
    mycloud.addEmitter('CO', 1.0e-4)

    lines = mycloud.lineLum('CO')
    
    nTrials += 1
    
    print nTrials
    
print time.time()-t0