''# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:38:08 2017

Generates contact map for each residue of the input trajectory file. Uses the alpha carbon as the location of each residue.

Input a PDB file or PDB trajectory file. Outputs average inter-residue distance.

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
def abDistance(atom1, atom2):
    return np.sqrt((atom1[0]-atom2[0])**2+(atom1[1]-atom2[1])**2+(atom1[2]-atom2[2])**2)       

#file =open("N_000_campari_traj.pdb", "r")
root = tkinter.Tk()
root.withdraw()
protein = input('What is the name of your protein?\n')

#Place all files requested into an array.
files = filedialog.askopenfilenames()

for file in files:
    #Data Arrays
    distanceMap = np.array([])
    model = []
    
    specificfile =open(file, "r")


    numOfResidues = 0
    m = 0
    
    #Determine the number of residues in the PDB file.
    for line in specificfile:
        
        if "ATOM" in line and "CA" in line:
            appendage = [m+1, int((str(line)[23:26]).strip()), str(line)[17:20], str(line)[13], returnMass(str(line)[13]), (str(line)[11:15]).strip(), float((str(line)[30:38]).strip()), float((str(line)[39:46]).strip()), float((str(line)[46:54]).strip())]
            model.append(appendage) 
    
        elif ("ENDMDL" in line) or ("END" in line):  
            #Make the distanceMap array.
            if m == 0:
                numOfResidues = model[-1][1]
                distanceMap = np.zeros((numOfResidues,numOfResidues))
            #Calculate distance of a residue from each residue in the chain.        
            for outer in range (0, numOfResidues):
                for inner in range (0, numOfResidues):
                    distanceMap[outer][inner] = distanceMap[outer][inner] + abDistance([model[inner][6],model[inner][7],model[inner][8]], [model[outer][6],model[outer][7],model[outer][8]])
            #Empty model array after data has been analyzed.
            model = []
            m = m + 1
            if m % 100 == 0:
                print ('Models processed: ' + str(m))
    distanceMap = distanceMap/m
    #If you want to limit the distance cutoff adjust cutOff variable (in angstroms)
    plt.imshow(distanceMap, cmap='gray', interpolation='nearest')
    plt.xlabel('Residue index')
    plt.ylabel('Residue index')
    plt.show()
    np.savetxt(str(file) + "_contactmap.txt", distanceMap, delimiter = ' ')            
    specificfile.close()

    


  
