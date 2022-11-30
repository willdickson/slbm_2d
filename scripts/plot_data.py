import sys
import numpy as np
import matplotlib.pyplot as plt

filename = sys.argv[1]

data = np.load(filename)

#data = data[:,:,::-1]
data = np.transpose(data,[0, 2, 1])
n = 4

p  = data[0,::n,::n]
ux = data[1,::n,::n]
uy = data[2,::n,::n]

fig, ax = plt.subplots(1,1)
ax.quiver(ux,uy)
#ax.axis('equal')
plt.show()


