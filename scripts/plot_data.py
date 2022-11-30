import sys
import numpy as np
import matplotlib.pyplot as plt

filename = sys.argv[1]

data = np.load(filename)

#data = data[:,:,::-1]
data = np.transpose(data,[0, 2, 1])
n = 1

p  = data[0,::n,::n]
ux = data[1,::n,::n]
uy = data[2,::n,::n]
print(uy[:,1:].shape, uy[:,:-1].shape)
w = (uy[1:, 1:] - uy[1: ,:-1]) - (ux[1:, 1:] - ux[:-1, 1:])

x = np.linspace(0,1,ux.shape[0])
y = np.linspace(0,1,uy.shape[0])
x, y = np.meshgrid(x,y)


fig, ax = plt.subplots(1,1)
w0 = np.absolute(w).max()
s = 0.03 
ax.pcolor(x,y,w, cmap='RdYlBu', vmin=-s*w0, vmax=s*w0)
#ax.quiver(ux,uy,scale=0.5)
ax.streamplot(x, y,ux,uy,color='k',linewidth=1.0, density=2.1)
ax.axis('equal')
plt.show()


