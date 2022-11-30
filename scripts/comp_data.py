import sys
import numpy as np
import matplotlib.pyplot as plt

colors = ['r', 'b', 'g', 'm', 'c']
#fig1, ax1 = plt.subplots(1,1)
#fig2, ax2 = plt.subplots(2,1)
fig3, ax3 = plt.subplots(1,1)

f_list = []
h_list = []

rho_list = []
ux_list = []
uy_list = []



for i, f in enumerate(sys.argv[1:]):

    print(f'filename: {f}')
    f_list.append(f)

    data = np.load(f)
    data = np.transpose(data,[0, 2, 1])
    
    rho  = data[0,:,:]
    ux = data[1,:,:]
    uy = data[2,:,:]

    rho_list.append(rho)
    ux_list.append(ux)
    uy_list.append(uy)

    print(ux[-2,:])

    
    #ax1.quiver(ux,uy, color=colors[i], scale=2.0)
    #ax2[i].pcolor(rho)
    h, = ax3.plot(ux[-2,:],'.-')
    #h_list.append(h)


dux = ux_list[0] - ux_list[1]
duy = uy_list[0] - uy_list[1]

print(dux[-2,:])

#ax3.quiver(dux, duy)




plt.show()


