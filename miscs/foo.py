import numpy as np

AvMax  = 30.0
dAvMin = 0.001
 
Av = 0.0
 
print Av
 
dAv = dAvMin
 
while Av <= AvMax + 1e-10:
     
    Av += dAv
      
    if Av < 1.0:
        dAv = 10.0**(int(np.log10(Av)))
         
    print Av