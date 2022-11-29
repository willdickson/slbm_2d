import sys
import numpy as np
import matplotlib.pyplot as plt

filename = sys.argv[1]

data = np.load(filename)

#data = data[:,:,::-1]
data = np.transpose(data,[0, 2, 1])

p  = data[0,:,:]
ux = data[1,:,:]
uy = data[2,:,:]

fig, ax = plt.subplots(1,1)
ax.quiver(ux,uy)
#ax.axis('equal')
plt.show()


