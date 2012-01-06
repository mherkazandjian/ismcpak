import h5py
from numpy import *
import numpy as np
Fx = []
nH = []
temp = []
new_Fx = []
new_nH = []
new_spe = []
new_tem = []
new_temp = []
filename = "deneme2.txt"
f = h5py.File('file1.hdf5', 'w')

count = 0
i = 1
for line in file(filename):        
    count = count+1
    if count>2:
        line = line.split()
        Fx.append(np.double(line[1]))
        nH.append(np.double(line[0]))
for i in range(len(Fx)):
    if float(nH[i]) == 8.00: 
        new_Fx.append(Fx[i])
    if nH[i] != nH[i-1]:
        new_nH.append(nH[i])
new_spe=["HI","HII","HeI","HeII","e","temp"]

#say = 0
#for line in file(filename):
#    say = say+1
#    if say>2:
#        line = line.split()
#        temp.append(line[6])
#        for k in range(len(temp)):
#            new_tem.append(temp[k])

dset = f.create_dataset('Parameter1', data=new_Fx)
dset = f.create_dataset('Parameter2', data=new_nH)

xdr = zeros((len(new_Fx), len(new_nH), len(new_spe)))
sayac = 0
for line in file(filename):
    sayac = sayac+1
    if sayac>2:
        line = line.split()
        for i in range(len(new_Fx)):
            if new_Fx[i]==np.double(line[1]):
                x=i
                break
        for j in range(len(new_nH)):
            if new_nH[j]==np.double(line[0]):
                y=j
                break
            
        xdr[x][y][0] = line[2]
        xdr[x][y][1] = line[3]
        xdr[x][y][2] = line[4]
        xdr[x][y][3] = line[5]
        xdr[x][y][4] = line[6]
        xdr[x][y][5] = line[7]
#        new_temp.append(line[7])

#dset = f.create_dataset('Temperature', data=new_temp, dtype='d')
dset = f.create_dataset('xdr', data=xdr, dtype ='d')
f.close()
#print dset
        
            
    
