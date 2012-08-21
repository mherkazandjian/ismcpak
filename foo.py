import numpy as np

Epsilon = 1e-10

AvMax  = 1.0
dAvMin = 0.01
 
AvArr = []

"""
i = 0; dAv = 0.0; Av = 0.0  # zeroth slab
AvArr.append(Av)
print i, Av
print '------------------'
"""

Av = 0.0

i = 1; dAv = dAvMin; Av += dAv  # first slab
AvArr.append(Av)
print i, Av
print '------------------'

while Av < 1.0 - Epsilon: # setting the Av up to Av = 1

    i += 1     
    Av += dAv
    AvArr.append(Av)
    print i, Av
    
    if abs(np.log10(Av)) % 1 <= 1e-14:
        dAv *= 10.0
        print 'increased dAv' 
    
dAv = 0.5
while Av < 10 - Epsilon :
    
    i += 1     
    Av += dAv
    AvArr.append(Av)
    print i, Av

dAv = 1.0
while Av < 30 - Epsilon :
    
    i += 1     
    Av += dAv
    AvArr.append(Av)
    print i, Av
    
    
print AvArr