# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 18:41:55 2020

Creates 2D histograms for Radius of gyration vs. Asphericity data.


@author: St. Peppers
"""
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm


#USER DEFINED VARIABLES
#Construct file
fontsize = 12
extension = '_campari_trajs_traj.pdb_contactmap.txt'
runs = ['N_000', 'N_001', 'N_002', 'N_003', 'N_004', 'N_005', 'N_006', 'N_007', 'N_008', 'N_009', 'N_010', 'N_011', 'N_012', 'N_013', 'N_014', 'N_015']
#Add the extension to each run name.
filenames = []
tempArray = [250,260,270,280,290,300,310,320,330,340,350,360,370,380,400,420]
tempInt = 0
#protein = input('What is the name of the construct?\n')

fig, ax = plt.subplots(4,4, sharex=True, sharey=True, figsize=(7,8))
fig.text(0.5, 0.08, r'Residue Index', ha='center', fontsize=15)
fig.text(0.05, 0.5, r'Residue Index', va='center', rotation='vertical', fontsize=15)

for i in range(0,len(runs)):
    data = np.loadtxt(runs[i]+extension)
    ax[int(i/4)][i%4].imshow(data, cmap='Spectral', interpolation='nearest')
    #ax[int(i/4)][i%4].set_ylabel(r'Asp $\delta$', fontsize = 15)
    #ax[int(i/4)][i%4].set_xlabel(r'R$_g$ ($\AA$)', fontsize = 15)
    ax[int(i/4)][i%4].set_title(str(tempArray[i])+'K')


