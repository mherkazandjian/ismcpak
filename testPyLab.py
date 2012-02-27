import matplotlib.pyplot as plt
import numpy as np

x = np.arange(6)
y = np.arange(6)

for i in xrange(10):
    if i==0:
        fig1 = plt.figure(figsize=(8, 6))
        ax1 = fig1.add_subplot(111)
        ax1.plot(x,y)
        fig1.show()
    else:
        ax1.plot(x+i,y+i)
        #p.set_data(z)

    print "step", i
    plt.pause(0.1)
