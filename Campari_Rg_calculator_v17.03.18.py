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
files = filedialog.askopenfilenames()


Macromolecule = input("What is the name of your macromolecule?\n")

tempArray = []
averageRgs = []
stdRgs = []

Temp = int(input("What is the lowest temperature RE?\n"))
Tincrement = input("What is the temperature increment?\n")

#Find the Rg for each model of each replica.  
for file in files:
    specificfile =open(file, "r")
    rgfile = open(str(Temp) + "_Rgs.txt", "w")
    processedmodels = 0
    model = []
    Rgs = []
    m = -1

    for line in specificfile:
        
        if "ATOM" in line:
            appendage = [m+1, int((str(line)[23:26]).strip()), str(line)[17:20], str(line)[13], returnMass(str(line)[13]), (str(line)[11:15]).strip(), float((str(line)[30:38]).strip()), float((str(line)[38:46]).strip()), float((str(line)[46:54]).strip())]
            model.append(appendage) 

        elif ("ENDMDL" in line) or ("END" in line):  
            
            #Determine Rg of last fully processed model. 
            xs = np.array(defX(model))
            ys = np.array(defY(model))
            zs = np.array(defZ(model))
            
            #Return the mass weighted center of the molecule.
 
            massArray = defMass(model)
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
            Rgs.append(newRg)
            rgfile.write("%s\n" % str(newRg))
            #End of Rg processing.
            if processedmodels % 100 == 0:
                print ("Models: ", processedmodels + 1) 
            #Empty model array after data has been analyzed.
            model = []
 
            processedmodels = processedmodels + 1
                        
    specificfile.close()
    rgfile.close()
    
    tempArray.append(Temp)
    averageRgs.append(np.average(Rgs))
    stdRgs.append(np.std(Rgs))
      
    
    plt.figure(figsize = (10,6))
    plt.hist(Rgs, 100, normed=1, facecolor='blue', alpha=0.75)
    plt.xlabel(r'$R_g$ ($\AA$)')
    plt.ylabel(r'P ($R_g$)')
    plt.title("R$_g$ distribution of " + Macromolecule + " at " + str(Temp) + "$^\circ$C")
    plt.savefig(str(Temp) + "_Rgs.png")
    plt.figure()
    plt.xlim(0,70)
    
    Temp = Temp + int(Tincrement)
  
plt.figure()
plt.scatter(tempArray,averageRgs)
plt.errorbar(tempArray,averageRgs,yerr=stdRgs, linestyle="None")
plt.xlabel(r'Temperature (K)')
plt.ylabel(r'$R_g$ ($\AA$)')
plt.savefig("RgvsTemp.png")