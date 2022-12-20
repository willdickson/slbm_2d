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
w = (uy[1:, 1:] - uy[1: ,:-1]) - (ux[1:, 1:] - ux[:-1, 1:])

x = np.linspace(0,ux.shape[1],ux.shape[1])
y = np.linspace(0,ux.shape[0],ux.shape[0])
x, y = np.meshgrid(x,y)


fig, ax = plt.subplots(1,1)
if 0:
    ax.pcolor(x,y,w, cmap='RdYlBu')
if 1:
    w0 = np.absolute(w).max()
    s = 0.5
    ax.pcolor(x,y,w, cmap='RdYlBu', vmin=-s*w0, vmax=s*w0)
    #ax.streamplot(x, y,ux,uy,color='k',linewidth=1.0, density=1.6)

    #ind = 100*np.arange(20)
    #ax.quiver(x[:,ind],y[:,ind],ux[:,ind],uy[:,ind])
    #print(ux[:,ind])

if 1:
    n = 2 
    tmp_ux = ux[::n,::n]
    tmp_uy = uy[::n,::n]
    tmp_x = x[::n,::n]
    tmp_y = y[::n,::n]
    ax.quiver(tmp_x,tmp_y,tmp_ux,tmp_uy,scale=20, headwidth=1, headlength=4, minlength=0)

if 0:
    for i in np.arange(20):
        s = 500
        y = np.arange(ux.shape[0])
        x0 = i*100
        k = ux.shape[1]//20
        print(k)
        x = x0 + s*ux[:,k*i]
        y = y[1:-1]
        x = x[1:-1]
    
        ax.plot([x0, x[0]], [y[0], y[0]], 'b')
        ax.plot([x0, x[0]], [y[-1],y[-1]], 'b')
        ax.plot(x,y,'b')
ax.axis('scaled')
plt.show()


