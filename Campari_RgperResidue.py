''# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:38:08 2017

This script formats Campari trajectorys (from a PDB format) into a custom format and calculates the Rg of each model.

@author: Erik Cook
"""

from tkinter import filedialog
import tkinter
import numpy as np

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
#Set the window size
window = 5
window = (window - 1)/2

#Find the Rg for each model of each replica.  

specificfile =open(file, "r")
processedmodels = 0
model = []
Rgs = []
m = -1
RgperRes = []
run = True
for line in specificfile:
    
    if "ATOM" in line:
        appendage = [m+1, int((str(line)[23:26]).strip()), str(line)[17:20], str(line)[13], returnMass(str(line)[13]), (str(line)[11:15]).strip(), float((str(line)[30:38]).strip()), float((str(line)[38:46]).strip()), float((str(line)[46:54]).strip())]
        model.append(appendage) 

    elif ("ENDMDL" in line) or ("END" in line):  
        resWindow = []
        if (run):
            lenOfProtein = model[-1][1]
            for i in range (0,lenOfProtein):
                RgperRes.append([])
            run = False
        else:
            None
        for residue in range (0,lenOfProtein):
            for atom in model:
                if (atom[1] >= (residue - window)) and (atom[1] <= (residue + window)):
                    resWindow.append(atom)
                else:
                    None      
            #Determine Rg of last fully processed model. 
            xs = np.array(defX(resWindow))
            ys = np.array(defY(resWindow))
            zs = np.array(defZ(resWindow))
            
            #Return the mass weighted center of the molecule.
     
            massArray = defMass(resWindow)
            cenx = np.sum(xs*massArray)/sum(massArray)
            ceny = np.sum(ys*massArray)/sum(massArray)
            cenz = np.sum(zs*massArray)/sum(massArray)
            
            xsquaredsum = 0
            ysquaredsum = 0
            zsquaredsum = 0
            
            i = 0
            while i < len(xs):
                xsquaredsum = xsquaredsum + (massArray[i]*(xs[i]-cenx)**2)
                i = i + 1
            
            i = 0  
            while i < len(ys):
                ysquaredsum = ysquaredsum + (massArray[i]*(ys[i]-ceny)**2)
                i = i + 1
                
            i = 0
            while i < len(zs):
                zsquaredsum = zsquaredsum + (massArray[i]*(zs[i]-cenz)**2)
                i = i + 1 
                
            i = 0
            newRg = np.sqrt((xsquaredsum+ysquaredsum+zsquaredsum)/(np.sum(massArray)))
            RgperRes[residue].append(newRg)
            #Empty the array.
            resWindow = []
        #End of Rg processing.
        if processedmodels % 100 == 0:
            print ("Models: ", processedmodels + 1) 
        else: 
            None
        #Empty model array after data has been analyzed.
        model = []
 
        processedmodels = processedmodels + 1
                    
specificfile.close()
Rgs = []
RgsStd = []
for i in RgperRes:
    Rgs.append(np.average(i))
    RgsStd.append(np.std(i))
plt.figure(figsize = (10,6))

Rgoutput = open('Rg_output280k.txt', 'w')
i = 0
while i < len(Rgs):
    Rgoutput.write(str(Rgs[i]) + ' ' + str(RgsStd[i]) + '\n')
    i = i + 1
Rgoutput.close()


plt.plot(np.arange(0,lenOfProtein), Rgs)
plt.xlabel(r'$R_g$ ($\AA$)')
plt.ylabel(r'P ($R_g$)')


plt.figure()
plt.xlim(0,70)
    

  
