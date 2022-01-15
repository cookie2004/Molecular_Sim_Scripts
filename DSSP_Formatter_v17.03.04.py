''# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:38:08 2017

Reads DSSP output from VMD and converts it into P(SS) vs. residue number.

@author: Erik Cook
"""


import tkinter
import numpy as np
import matplotlib.pyplot as plt


     
root = tkinter.Tk()
root.withdraw()

file = "N_007_SS.txt"
#Place all files requested into an array.
fname = open("N_007_SS.txt", 'r')
lines=fname.readlines()
residues = int((np.array(lines[6].split()))[2])
frames = int((np.array(lines[5].split()))[2])
SS_container = []

Helix = []
Strand = []
Coil = []



residue = 0
residueArray = []
# The following variable represents the number of header lines in the incoming file. Adjust accordingly.
offset = 9
residue = 1
i = 10
frame = 0
##This variable 'frames' is here for testing. Remove when you want to add the real thing.


while residue <= residues:
    print (residue)
    residueArray.append(residue)
    frame = 0
    
    while frame < frames:
        SS_container.append(lines[offset + residue + frame*residues-1].split()[4])
        frame = frame + 1
    H = 0
    C = 0
    S = 0  
    for s in SS_container:
        if ('H' in s) or ('G' in s) or ('I' in s):
            H = H + 1
        elif ('S' in s) or ('T' in s) or ('C' in s):
            C = C + 1
        elif ('E' in s) or ('B' in s):
            S = S + 1
        
    Helix.append(H/frames)
    Strand.append(S/frames) 
    Coil.append(C/frames)   
    residue = residue + 1
    SS_container = []   
    
#for line in file:
#    #Determine the number of residues.
#    
#    for x in line:
        
fname.close()
nStrand = []
for i in Strand:
    i = i * -1
    nStrand.append(i)
    


figure = plt.bar(residueArray, Helix, color = "black", edgecolor = "white", width = 1)

plt.title("PDX1NTAD at 300K")

figure2 = plt.bar(residueArray, nStrand, color = "red")
plt.ylim(-1,1)
plt.xlim(1,residues)
plt.xlabel("Residue index")
plt.ylabel("$P$")
s
plt.figure()
figure3 = plt.bar(residueArray, Coil)
plt.ylim(0,1)
plt.title("Coil content")
plt.figure()