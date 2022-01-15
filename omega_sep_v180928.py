''# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:38:08 2017

This script formats Campari trajectorys (from a PDB format) into a custom format and calculates the cis/trans residue abundance.


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
    
def dihedral(p):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = np.array(p[0])
    p1 = np.array(p[1])
    p2 = np.array(p[2])
    p3 = np.array(p[3])

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))
 

#Place all files requested into an array.
file = open('TEST.pdb')
print('Models will be split into trans-P6 and cis-P6.')
#Make cis and trans subsets.
Macromolecule = 'CTD'
cis_models = open(Macromolecule + "cis_models.pdb", 'w')
trans_models = open(Macromolecule + "trans_models.pdb", 'w')

tempArray = []
averageRgs = []
stdRgs = []
residue = 21

processedmodels = 0
model = []
model_str = []
m = -1

cisocc = []
for i in range(0,residue):
    cisocc.append(0)

model = []
for i in range(0,residue):
    model.append([0,0,0,0])
    
for line in file:
    
    if ' C ' in line or ' O ' in line or ' N ' in line or ' HN' in line or (' PRO ' in line and ' CD ' in line):
        appendage = [m+1, int((str(line)[23:26]).strip()), str(line)[17:20], str(line)[13], returnMass(str(line)[13]), (str(line)[11:15]).strip(), float((str(line)[30:38]).strip()), float((str(line)[38:46]).strip()), float((str(line)[46:54]).strip())]
        

        if appendage[5] == 'O':
            model[appendage[1]-1][0] = [appendage[6], appendage[7], appendage[8]]
        elif appendage[5] == 'C':
            model[appendage[1]-1][1] = [appendage[6], appendage[7], appendage[8]]
        elif appendage[5] == 'N':
            model[appendage[1]-1][2] = [appendage[6], appendage[7], appendage[8]]
        elif appendage[5] == 'HN':
            model[appendage[1]-1][3] = [appendage[6], appendage[7], appendage[8]]
        elif appendage[5] == 'CD':
            model[appendage[1]-1][3] = [appendage[6], appendage[7], appendage[8]]      
        else:
            None
    if "ATOM" in line:
        model_str.append(line)
    elif ("ENDMDL" in line) or ("END" in line):  
        
        #Calculate the omega bond angles. From -180
        for i in range(0,residue-1):
            p = model[i][0],model[i][1], model[i+1][2], model[i+1][3]
            if dihedral(p) <90 and dihedral(p) > -90:
                cisocc[i] += 1
                if int(i) == 11: 
                    for el in model_str:
                        cis_models.write(el)
                    cis_models.write('ENDMDL\n')
                else:
                    for el in model_str:
                        trans_models.write(el)
                    trans_models.write('ENDMDL\n')
            else: 
                None
        
        #Empty model array to prepare for the next iteration. 
        model = []
        model_str = []
        for i in range(0,residue):
            model.append([0,0,0,0])
        
        if appendage[0]%100 == 0:
            print ("Models processed: " + str(appendage[0]))
        else:
            None
        m +=1
    else:
        None
        
        
#Convert occurances of cis residues into probabilities.

for i in range(0,len(cisocc)):
    cisocc[i] = (cisocc[i]/(m+1))*100
    
cis_models.close()
trans_models.close()

#Construct plots.

dataX = np.arange(1,len(cisocc)+1)
plt.figure(figsize=(7,2))
plt.xlim(0,residue+1)
plt.xlabel('Residue Index')
plt.ylabel('% cis')
plt.bar(dataX+1, cisocc, color = 'black')
plt.savefig('%cis.pdf')   

