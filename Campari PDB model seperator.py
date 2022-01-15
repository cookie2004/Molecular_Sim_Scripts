''# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:38:08 2017

This script formats Campari trajectorys (from a PDB format) into a custom format and calculates the Rg of each model.

@author: Erik Cook
"""

from tkinter import filedialog
import tkinter
       
#file =open("N_000_campari_traj.pdb", "r")
root = tkinter.Tk()
root.withdraw()

#Place all files requested into an array.
files = filedialog.askopenfilenames()

#Find the Rg for each model of each replica.  



model = []

m = 0
for file in files: 
    specificfile =open(file, "r")
    for line in specificfile:
        
        if "ATOM" in line and "CA" in line:
            model.append(line)
    
        elif ("ENDMDL" in line) or ("END" in line):  
            
            output = open(str(m)+'eom.pdb', 'w')
            for i in model:
                output.write(i)
            output.close()
            model = []
            m = m + 1
            if m % 100 == 0:
                print ('Models completed: ' + str(m))
    specificfile.close()                       
