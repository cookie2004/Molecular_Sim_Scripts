''# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:38:08 2017

This script formats Campari trajectorys (from a PDB format) into a custom format composing of only the backbone atoms.

@author: Erik Cook
"""

from tkinter import filedialog
import tkinter


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

Otherfiles = ['N_001_campari_traj.pdb','N_002_campari_traj.pdb','N_003_campari_traj.pdb','N_004_campari_traj.pdb','N_005_campari_traj.pdb','N_006_campari_traj.pdb','N_007_campari_traj.pdb','N_008_campari_traj.pdb','N_009_campari_traj.pdb','N_010_campari_traj.pdb','N_011_campari_traj.pdb','N_012_campari_traj.pdb','N_013_campari_traj.pdb','N_014_campari_traj.pdb','N_015_campari_traj.pdb','N_016_campari_traj.pdb','N_017_campari_traj.pdb','N_018_campari_traj.pdb','N_019_campari_traj.pdb']

#Place all files requested into an array.
files = ['N_007_campari_traj.pdb']

fileName = input("Name your output file.\n")
fileNumber = 0
processedmodels = 0
#Find the Rg for each model of each replica.  
for file in files:
    pdbTraj =open(file, "r")
    Cafile = open(str(fileName) + "_" + str(fileNumber) + "_Ca.txt", "w")
    
    m = 1
    
    for line in pdbTraj:
        if ("ATOM" in line) and (" N " in line or "CA" in line or " C " in line):
            
            Cafile.write(str(m) + " " + str(line)[23:26].strip() + " " + str(line)[17:20] + " " + str(line)[13] + " " + str(returnMass(str(line)[13])) + " " + str(line)[11:15].strip() + " " + str(line)[30:38].strip() + " " + str(line)[39:46].strip() + " " + str(line)[46:54].strip() + "\n")
    
        elif ("ENDMDL" in line) or ("END" in line):  
            if processedmodels % 100 == 0:
                print ("Processed models: ", processedmodels + 1)  
            
            Cafile.write("END\n")
            m = m + 1 

            processedmodels = processedmodels + 1
                        
    pdbTraj.close()
    Cafile.close()
    fileNumber = fileNumber + 1


      
    
