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
from numpy import linalg as LA
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
    asfile = open(str(Temp)+ "_Asphericity.txt", "w")
    processedmodels = 0
    model = []
    Rgs = []
    asphericity = []
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

            gTen = np.empty([3,3])
            for atom in range(0,len(xs)):
                Sxx = (xs[atom]-cenx)**2
                Sxy = (xs[atom]-cenx)*(ys[atom]-ceny)
                Syy = (ys[atom]-ceny)**2
                Syz = (ys[atom]-ceny)*(zs[atom]-cenz)
                Szz = (zs[atom]-cenz)**2
                Sxz = (xs[atom]-cenx)*(zs[atom]-cenz)
                gTen = gTen + np.array([[Sxx, Sxy, Sxz], [Sxy, Syy, Syz], [Sxz, Syz, Szz]])
                
            gTen = gTen/len(xs)
            
            #Diagonalize the gyration tensor
            pC = LA.eig(gTen)
            pM = np.sort(pC[0])
            #Determine asphericity and append it to the array.
            newAsp = 1-3*((pM[0]*pM[1]+pM[1]*pM[2]+pM[2]*pM[0])/(np.sum(pM)*np.sum(pM)))
            

            newRg = np.sqrt(np.sum(pM))
            #Check new asphericity value to see if it is valid. If not, do not append it to the asphericity data set.
            if (newAsp <= 1 and newAsp >= 0):
                Rgs.append(newRg)
                asphericity.append(newAsp)
            else:
                None
            rgfile.write("%s\n" % str(newRg))
            asfile.write("%s\n" % str(newAsp))
            
            #End of Rg processing.
            if processedmodels % 100 == 0:
                print ("Models processed: ", processedmodels + 1) 
           
            #Empty model array after data has been analyzed.
            model = []
            processedmodels = processedmodels + 1
            gTen = np.empty([3,3])
        else:
            None
        
        
    specificfile.close()
    rgfile.close()
    asfile.close()
    tempArray.append(Temp)
    averageRgs.append(np.average(Rgs))
    stdRgs.append(np.std(Rgs))
      
    
    plt.figure(figsize = (6,6))
    plt.hist(Rgs, 50, normed=1, facecolor='blue', alpha=0.75)
    plt.xlabel(r'$R_g$ ($\AA$)',fontsize = 15)
    plt.ylabel(r'P ($R_g$)',fontsize = 15)
    plt.title("R$_g$ distribution of " + Macromolecule + " at " + str(Temp) + "$^\circ$C")
    plt.savefig(str(Temp) + "_Rgs.pdf")
 
    plt.figure(figsize=(6,6))
    x = np.arange(0,len(asphericity))
    plt.title("Asphericity distribution of " + Macromolecule + " at " + str(Temp) + "$^\circ$C")
    plt.scatter(Rgs, asphericity, c=x)
    plt.ylabel('$\delta$',fontsize = 15)
    plt.xlabel(r'$R_g$ ($\AA$)',fontsize = 15)
    plt.ylim(0,1)
    plt.tight_layout()
    plt.savefig(str(Temp) + "_Asphericity.pdf")
    
    Temp = Temp + int(Tincrement)
  
plt.figure(figsize=(6,6))
plt.scatter(tempArray,averageRgs)
plt.errorbar(tempArray,averageRgs,yerr=stdRgs, linestyle="None")
plt.xlabel(r'Temperature (K)')
plt.ylabel(r'$R_g$ ($\AA$)')
plt.savefig("RgvsTemp.pdf")



