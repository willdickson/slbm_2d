import sys
import numpy as np
import matplotlib.pyplot as plt

mesh_id = np.load(sys.argv[1])

fig, ax = plt.subplots(1,1)

ax.pcolor(mesh_id)
plt.show()
