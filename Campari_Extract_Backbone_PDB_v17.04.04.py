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
file = 'N_007_campari_traj.pdb'

fileName = input("Name your output file.\n")
model = 0
#Find the Rg for each model of each replica.  

pdbTraj =open(file, "r")
Cafile = open(str(fileName) + "_Ca.pdb", "w")

for line in pdbTraj:
    if "CRYST" in line:
        if model % 100 == 0:
            Cafile.write(line)
    if "TITLE" in line:
        if model % 100 == 0:
            Cafile.write(line)
    if "ATOM" in line and ("CA" in line or "N" in line or " C " in line):
        if model % 100 == 0:
            Cafile.write(line)

    elif ("ENDMDL" in line) or ("END" in line):  
         
        if model % 100 == 0:
            Cafile.write(line)
        model = model + 1
        if model % 100 == 0:
            print (model)

                    
pdbTraj.close()
Cafile.close()



      
    
