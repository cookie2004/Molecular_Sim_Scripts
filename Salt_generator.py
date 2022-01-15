# -*- coding: utf-8 -*-
"""
Created on Thu May  4 11:12:48 2017

@author: Erik Cook
"""

import numpy as np

radius = input('What is the radius of your spherical droplet (in angstoms)?\n')
V = 4/3*np.pi*float(radius)**3
concentration = input('Input your concentration of KCl needed.\n')
numOfAtoms = float(concentration) * (1*10**(-27))*V*(6.022*(10**23))

file = open(str(concentration)+'M_KCl.txt','w')

x = 0
while x < numOfAtoms:  
    file.write('K+\n')
    x = x + 1
x = 0
while x < numOfAtoms:
    file.write('CL-\n')
    x = x + 1
file.close()