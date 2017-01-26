
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 16:17:36 2017

@author: apm13
plot heat map
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import griddata

filename=sys.argv[-1]
data=np.genfromtxt(filename, delimiter=' ')
X_dat=data[:,0]
Y_dat=data[:,1]
E_dat=data[:,6]
V_dat=data[:,7]

# Convert from pandas dataframes to numpy arrays
X, Y, E, V, = np.array([]), np.array([]), np.array([]), np.array([])
for i in range(len(X_dat)):
        X = np.append(X,X_dat[i])
        Y = np.append(Y,Y_dat[i])
        E = np.append(E,E_dat[i])
        V = np.append(V,V_dat[i])

# create x-y points to be used in heatmap
xi = np.linspace(X.min(),X.max(),100)
yi = np.linspace(Y.min(),Y.max(),100)

# E,V are matrices of x-y values
Ei = griddata((X, Y), E, (xi[None,:], yi[:,None]), method='cubic')
Vi = griddata((X, Y), V, (xi[None,:], yi[:,None]), method='cubic')

plt.figure(1)
CS = plt.contourf(xi, yi, Ei, 10, cmap=plt.cm.rainbow)
plt.savefig("Eheat.pdf")
plt.colorbar()
plt.show()
plt.close()

plt.figure(2)
CS = plt.contourf (xi, yi, Vi, 10, cmap=plt.cm.rainbow)
plt.savefig("Vheat.pdf")
plt.colorbar()
plt.show()
plt.close()