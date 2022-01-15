''# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:38:08 2017

This script formats Campari trajectorys (from a PDB format) into a custom format and calculates the cis/trans residue abundance.


@author: Erik Cook
"""

from tkinter import filedialog
import tkinter
import numpy as np
import sys
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
 
fontsize = 12
extension = '_campari_trajs_traj.pdb'
runs = ['N_000', 'N_001', 'N_002', 'N_003', 'N_004', 'N_005', 'N_006', 'N_007', 'N_008', 'N_009', 'N_010', 'N_011', 'N_012', 'N_013', 'N_014', 'N_015']
#Add the extension to each run name.
filenames = []
tempArray = [250,260,270,280,290,300,310,320,330,340,350,360,370,380,400,420]
tempInt = 0
#protein = input('What is the name of the construct?\n')

fig, ax = plt.subplots(len(tempArray),1, sharex=True)
fig.set_figheight(len(tempArray)*1.5)
fig.subplots_adjust(hspace=0.15)

fig.text(0.5, 0.08, r'Residue Index', ha='center', fontsize=15)
fig.text(0.08, 0.5, r'P(Cis-Pro)', va='center', rotation='vertical', fontsize=15)

for run in runs:
    filename = run+extension
    #Place all files requested into an array.
    file = open(filename)
    residues = 112
    sys.stdout.write("Processing run " + str(run)+'\n')
    #Macromolecule = input("What is the name of your macromolecule?\n")
    
    tempArray = []
    averageRgs = []
    stdRgs = []
    
    
    #Find the Rg for each model of each replica.  
    processedmodels = 0
    model = []
    m = -1
    
    cisocc = []
    for i in range(0,residues):
        cisocc.append(0)
    
    model = []
    for i in range(0,residues):
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
                
        elif ("ENDMDL" in line) or ("END" in line):  
            
            #Calculate the omega bond angles. From -180
            for i in range(0,residues-1):
                p = model[i][0],model[i][1], model[i+1][2], model[i+1][3]
                if dihedral(p) <90 and dihedral(p) > -90:
                    cisocc[i] += 1
                else: 
                    None
            
            #Empty model array to prepare for the next iteration. 
            model = []
            for i in range(0,residues):
                model.append([0,0,0,0])
            
            if appendage[0]%100 == 0:
                sys.stdout.write("\rModels processed: " + str(appendage[0]))
                sys.stdout.flush()
            else:
                None
            m +=1
        else:
            None
    sys.stdout.write("\nDone!\n")     
    #Convert occurances of cis residues into probabilities.
    
    for i in range(0,len(cisocc)):
        cisocc[i] = (cisocc[i]/(m+1))*100
        
    np.savetxt(run+'cis_pro.txt', cisocc)
    #Construct plots.
    
    dataX = np.arange(0,len(cisocc))
    ax[tempInt].bar(dataX, cisocc, color = 'black')
    ax[tempInt].set_ylim(0,50)
    tempInt += 1


