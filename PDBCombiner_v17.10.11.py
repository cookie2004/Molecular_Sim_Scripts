# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 13:59:46 2017

This scripts combines PDB files into one file.
@author: Erik Cook
"""

''# -*- coding: utf-8 -*-

from tkinter import filedialog
import tkinter   
#file =open("N_000_campari_traj.pdb", "r")
root = tkinter.Tk()
root.withdraw()

#Place all files requested into an array.
files = filedialog.askopenfilenames()

#Find the Rg for each model of each replica.  

model = []

output = input('Name your output file.\n')
outFile = open(output + '.pdb', "w")

m = 0
for file in files: 
    outFile.write("MODEL" + str(m)+"\n")
    specificfile =open(file, "r")
    for line in specificfile:
      line = line.replace('\n','')
      line = line.strip()
      if ("END" in line) or not(line) or len(line)==0 or not(line.rstrip()):
          continue
      else:
          outFile.write(line + "\n")
    if m % 100 == 0:
        print ('Models completed: ' + str(m))
    else:
        None
    outFile.write('ENDMDL\n')
    specificfile.close()      
    m = m + 1


outFile.close()                 
