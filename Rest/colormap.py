#%%
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
size = 256
size = range(35,size,1)

cmap = cm.nipy_spectral(range(35,256,1))
cmap = np.interp(cmap, 256, 256)
size = np.interp(size, 256, 256)
print(cmap)
plt.scatter(size, size, c=cmap)
plt.show()