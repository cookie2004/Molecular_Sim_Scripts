''# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:38:08 2017

This script formats Campari trajectorys (from a PDB format) into a custom format and calculates the Rg of each model.

@author: Erik Cook
"""

from tkinter import filedialog
import tkinter
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from PIL import Image
from io import BytesIO
 
def defX(matrix):
    return [row[6] for row in matrix]

def defY(matrix):
    return [row[7] for row in matrix]

def defZ(matrix):
    return [row[8] for row in matrix]

def defMass(matrix):
    return [row[4] for row in matrix]

def returnMass(atom):
    if "C" in atom:
        return 12.0107
    if "N" in atom:
        return 14.0067
    if "O" in atom:
        return 15.9994
    if "H" in atom:
        return 1.00794
    if "S" in atom:
        return 32.065
        
#file =open("N_000_campari_traj.pdb", "r")
root = tkinter.Tk()
root.withdraw()


#Place all files requested into an array.
file = filedialog.askopenfilename()

tempArray = []
averageRgs = []
stdRgs = []


#Find the Rg for each model of each replica.  

specificfile =open(file, "r")
processedmodels = 0
model = []

Rgs = []

for i in range(0,82):
    Rgs.append([0])

for j in range(0,82):
    for i in range(0,81):
        Rgs[j].append([0,0])

        
        
m = -1

run = True
for line in specificfile:
    
    if "ATOM" in line:
        appendage = [m+1, int((str(line)[23:26]).strip()), str(line)[17:20], str(line)[13], returnMass(str(line)[13]), (str(line)[11:15]).strip(), float((str(line)[30:38]).strip()), float((str(line)[38:46]).strip()), float((str(line)[46:54]).strip())]
        model.append(appendage) 
        
        
    elif ("ENDMDL" in line) or ("END" in line): 
        filteredmodel = []
        for i in model:
            if 'CA' in i:
                filteredmodel.append(i)
            else: 
                None
        ij = []
        for i in filteredmodel:
            for j in filteredmodel:
                

                r = math.sqrt(((i[6] - j[6])**2 + (i[7] - j[7])**2 + (i[8]-i[8])**2))
                #r = r/(i[1]/j[1])
                Rgs[i[1]][j[1]][1] += r
                Rgs[i[1]][j[1]][0] += 1
                
                
                #Empty the array.
                
        #End of Rg processing.
        
        if processedmodels % 100 == 0:
            print ("Models: ", processedmodels + 1) 
        else: 
            None
        #Empty model array after data has been analyzed.
        model = []
 
        processedmodels = processedmodels + 1
                    
#Process the array
Rgs.pop(0)
for i in Rgs:
    i.remove(0)

for i in Rgs:
    for j in i:
        j[1] = j[1]/j[0]
        j[0] = j[0]/j[0]

#Configure X-axis       
for i in Rgs:
    x = 1 
    for j in i:
        j[0] = x
        x += 1
x = 0
while x < len(Rgs):
    for i in Rgs[x]:
        i[0] = abs(i[0]-(x+1))
    x += 1
x= 203    
for i in Rgs:
    plt.figure(figsize=(4,4))
    plt.xlim(0,83)
    data = np.array(i)
    plt.scatter(data[:,0], data[:,1], color = 'black')
    plt.title('Residue: ' + str(x))
    plt.xlabel('|i-j|')
    plt.ylabel(r'<R$_{ij}$>')
    x += 1
