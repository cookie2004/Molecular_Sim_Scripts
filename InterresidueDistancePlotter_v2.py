''# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:38:08 2017

This script formats Campari trajectorys (from a PDB format) into a custom format composing of only the backbone atoms.

@author: Erik Cook
"""

from tkinter import filedialog
import tkinter
import numpy as np
import matplotlib.pyplot as plt

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
#correlation time for the electron nuclear interaction
tau = 0.0000000207
#Larmor frequency of the nuclear spin (proton, carbon, etc.) For proton, omega = 500,000,000Hz. For carbon, omega = 125,700,000Hz. 
H1omega = 498.186000
C13omega = 125.7
K = 1.23e-32        
t = 0.009
H1R2 = 0.00001
C13R2 = 0.0001

def R2SP(tau, omega, K, r):
    r = r*10e-8
    b = 4*tau+((3*tau)/(1+omega**2*tau**2))
    return (K * b)/(r**6)

def Iratio(R2, R2S, t):
    return (R2**(-R2S*t))/(R2 + R2S)
    
def distance(K, tau, omega, R2S):
    a = K/R2S
    b = 4*tau+((3*tau)/(1+omega**2*tau**2))
    return (a*b)**(1/6)
    
#file =open("N_000_campari_traj.pdb", "r")
root = tkinter.Tk()
root.withdraw()


#Place all files requested into an array.
files = filedialog.askopenfilenames()

fileName = input("Name your output file.\n")
fileNumber = 0
processedmodels = 0


#Find the Rg for each model of each replica.  
for file in files:
    pdbTraj =open(file, "r")
    model = []
    m = 1
    Cys = []
    distanceArray = np.array([])
    aveDistance = np.array([])
    for line in pdbTraj:
        
        if "ATOM" in line:
            appendage = [float((str(line)[30:38]).strip()), float((str(line)[39:46]).strip()), float((str(line)[46:54]).strip())]
            model.append(appendage)
        
            #Use the location of the sulpher group as the location of the MTSL label.
        if "CYS" in line and "CA" in line:
            Cys = [float((str(line)[30:38]).strip()), float((str(line)[39:46]).strip()), float((str(line)[46:54]).strip())]
        
        elif ("ENDMDL" in line) or ("END" in line) or not(line):  
            
            for i in model:
                dis = np.sqrt(((i[0]-Cys[0])**2)+((i[1]-Cys[1])**2)+((i[2]-Cys[2])**2))
                distanceArray = np.append(distanceArray, dis)
 
            if len(aveDistance) == 0:
                aveDistance = distanceArray
            else:
                aveDistance = aveDistance + distanceArray
            
            #Empty distanceArray and model array.
            distanceArray = np.array([])
            model = []

            if processedmodels % 100 == 0:
                print (processedmodels)    
            processedmodels = processedmodels + 1

    aveDistance = aveDistance/processedmodels
    np.savetxt(str(fileName) + '.txt', aveDistance)
    residue = np.arange(1, len(aveDistance)+1, 1)     
    pdbTraj.close()
    
    C13R2SP = np.array([])
    H1R2SP = np.array([])
    
 
    H1R2SP = R2SP(tau, H1omega, K, aveDistance)
    H1ratio = Iratio(H1R2, H1R2SP, t)
    H1ratio = H1ratio - np.min(H1ratio)
    H1ratio = H1ratio/np.max(H1ratio)
    
    C13R2SP = R2SP(tau, C13omega, K, aveDistance)
    C13ratio = Iratio(C13R2, C13R2SP, t)
    C13ratio = C13ratio - np.min(C13ratio)
    C13ratio = C13ratio/np.max(C13ratio)
    
    plt.plot(residue, C13R2SP)    
    fig,ax=plt.subplots(figsize = (10,6))
    figure = plt.plot(residue, aveDistance)
    plt.title("PDX1C at 350K", fontsize = "20")
    fileNumber = fileNumber + 1    



  
    