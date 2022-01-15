# -*- coding: utf-8 -*-
"""
Created on Sun May  3 11:51:04 2020

@author: St. Peppers
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import savgol_filter
from sklearn.decomposition import PCA as sk_pca
from sklearn.preprocessing import StandardScaler
from sklearn import svm
from sklearn.cluster import KMeans
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
from sklearn.datasets.samples_generator import make_blobs
import scipy.stats as st
import scipy.constants as ct

def BoltzmannDist(E, T, S):
    return 2 * np.sqrt(E/ct.pi)*(1/(ct.Boltzmann*T))**(S) * np.exp(-E/(ct.Boltzmann*T))

def normalizeArray(array):
    output = array-(np.min(array)*0.999)
    output = output/np.max(output)
    return output

def redefine(array, define):
    output=[]
    for i in array:
        i_array = []
        for j in i:
            if j < define:
                i_array.append(j)
            else:
                i_array.append(define)
        output.append(i_array)
    return np.array(output)
temps = [250.0,260.0,270.0,280.0,290.0,300.0,310.0,320.0,330.0,340.0,350.0,360.0,370.0,380.0,400.0,420.0]
folders_wtRD = ["190925_RDwt_2","190926_wtRD"]
folders_RDmut = ["190925_RD_E415K_E418K_2","190926_RD_E415K_E418K"]



# for folder in folders_wtRD:
#     for i in range (5,6):
#     #Perform WHAM construction. 



data=[[],[]]
for i in range (0,16):
    data[0].append(np.loadtxt(str(i)+str(folders_wtRD[0]+'_PCs.txt'))[:,0])
    data[1].append(np.loadtxt(str(i)+str(folders_wtRD[0]+'_PCs.txt'))[:,1])
    
data1 = np.loadtxt('5190925_RDwt_2_PCs.txt')
data2 = np.loadtxt('5190926_wtRD_PCs.txt')
RDmut1 = np.loadtxt('5190926_RD_E415K_E418K_PCs.txt')
RDmut2 = np.loadtxt('5190925_RD_E415K_E418K_2_PCs.txt')

wtRD_PC1 = np.append(data1[:,0], data2[:,0],axis=0)
wtRD_PC2 = np.append(data1[:,1], data2[:,1], axis=0)

RDmut_PC1 = np.append(RDmut1[:,0], RDmut2[:,0], axis=0)
RDmut_PC2 = np.append(RDmut1[:,1], RDmut2[:,1], axis=0)

plt.figure(figsize=(5,5))
bins = plt.hist2d(wtRD_PC1, wtRD_PC2,500)

# plt.figure(figsize=(5,5))
# PC1_BD = []
# for i in data[0]:
#     PC1_BD.append(plt.hist(i, bins[1]))


# plt.figure(figsize=(5,5))
# wtRDbins = plt.hist2d(wtRD_PC1, wtRD_PC2,500, cmap=plt.cm.Spectral, norm=LogNorm())
# plt.xlim(-25,25)
# plt.ylim(-25,25)

# plt.figure(figsize=(5,5))
# plt.scatter(wtRD_PC1,wtRD_PC2)


#3D density plot
n_components = 3
X, truth = make_blobs(n_samples=300, centers=n_components, 
                      cluster_std = [2, 1.5, 1], 
                      random_state=42)

#Extract x and y
x = wtRD_PC1
y = wtRD_PC2
deltaX = (max(x) - min(x))/10
deltaY = (max(y) - min(y))/10
xmin = min(x) - deltaX
xmax = max(x) + deltaX
ymin = min(y) - deltaY
ymax = max(y) + deltaY
print(xmin, xmax, ymin, ymax)
# Create meshgrid
xx, yy = np.mgrid[-25:10:100j, -15:15:100j]

plt.figure(figsize=(7, 7))
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([x, y])
kernel = st.gaussian_kde(values)
kernel.set_bandwidth(bw_method=kernel.factor / 12.)
f = normalizeArray(np.reshape(kernel(positions).T, xx.shape))
f_log = -np.log(f)
f_log = redefine(f_log, 10)
fig = plt.figure(figsize=(8,8))
ax = fig.gca()
ax.set_xlim(-25,10)
ax.set_ylim(-15,15)
cfset = ax.contourf(xx, yy, f_log, cmap='Spectral')
ax.imshow(np.rot90(f), cmap='Spectral', extent=[xmin, xmax, ymin, ymax])
cset = ax.contour(xx, yy, f_log, colors='k')
#ax.clabel(cset, inline=0, fontsize=10)
ax.set_xlabel('PC1', fontsize=15)
ax.set_ylabel('PC2', fontsize=15)
plt.title('wtRD Folding landscape', fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.annotate('X - State 1', (-10,-5), fontsize=20)
plt.annotate('X - State 2', (-5,2.5), fontsize=20)
plt.annotate('X - State 3', (-12,6), fontsize=20)
plt.savefig('wtRD_FL.pdf', dpi=300)

fig = plt.figure(figsize=(13, 7))
ax = plt.axes(projection='3d')
surf = ax.plot_surface(xx, yy, f_log, rstride=1, cstride=1, cmap='Spectral', edgecolor='none')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PMF (kcal/mol)')
ax.set_title('wtRD Folding landscape')
fig.colorbar(surf, shrink=0.5, aspect=5) # add color bar indicating the PDF
ax.view_init(20, 15)
# Get rid of the ticks
ax.set_xticks([]) 
ax.set_yticks([]) 
ax.set_zticks([])
# Extract x and y
x = RDmut_PC1
y = RDmut_PC2
deltaX = (max(x) - min(x))/10
deltaY = (max(y) - min(y))/10
xmin = min(x) - deltaX
xmax = max(x) + deltaX
ymin = min(y) - deltaY
ymax = max(y) + deltaY
print(xmin, xmax, ymin, ymax)
# Create meshgrid
xx, yy = np.mgrid[-20:0:100j, -15:15:100j]


#Plot the RD mutant data.
#Extract x and y
x = RDmut_PC1
y = RDmut_PC2
deltaX = (max(x) - min(x))/10
deltaY = (max(y) - min(y))/10
xmin = min(x) - deltaX
xmax = max(x) + deltaX
ymin = min(y) - deltaY
ymax = max(y) + deltaY
print(xmin, xmax, ymin, ymax)
# Create meshgrid
xx, yy = np.mgrid[-25:10:100j, -15:15:100j]

plt.figure(figsize=(7, 7))
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([x, y])
kernel = st.gaussian_kde(values)
kernel.set_bandwidth(bw_method=kernel.factor / 12.)
f2 = normalizeArray(np.reshape(kernel(positions).T, xx.shape))
f2_log = -np.log(f2)
f2_log = redefine(f2_log, 10)
fig = plt.figure(figsize=(8,8))
ax = fig.gca()
ax.set_xlim(-25,10)
ax.set_ylim(-15,15)
cfset = ax.contourf(xx, yy, f2_log, cmap='Spectral')
ax.imshow(np.rot90(f), cmap='Spectral', extent=[xmin, xmax, ymin, ymax])
cset = ax.contour(xx, yy, f2_log, colors='k')
#ax.clabel(cset, inline=1, fontsize=10)
ax.set_xlabel('PC1', fontsize=15)
ax.set_ylabel('PC2', fontsize=15)
plt.title('RD E415K, E418K Folding landscape', fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.annotate('X - State 1', (-10,-5), fontsize=20)
plt.annotate('X - State 2', (-5,2.5), fontsize=20)
plt.annotate('X - State 3', (-12,6), fontsize=20)
plt.savefig('RDmut_FL.pdf', dpi=300)

fig = plt.figure(figsize=(13, 7))
ax = plt.axes(projection='3d')
surf = ax.plot_surface(xx, yy, f2_log, rstride=1, cstride=1, cmap='Spectral', edgecolor='none')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PMF kcal/mol')
ax.set_title('RD E415K, E418K Folding landscape')
fig.colorbar(surf, shrink=0.5, aspect=5) # add color bar indicating the PDF
ax.view_init(20, 15)
# Get rid of the ticks
ax.set_xticks([]) 
ax.set_yticks([]) 
ax.set_zticks([])
plt.savefig('colorbar.pdf', dpi=300)


