# -*- coding: utf-8 -*-
"""
Created on Mon May  1 15:40:32 2017

Determines the number of Na+ or Cl- ions required to balance net charge of system.

@author: Erik Cook
"""

from tkinter import filedialog
import tkinter
import numpy as np
import scipy as sp

root = tkinter.Tk()
root.withdraw()

#Import .key file containing protein sequence.

file = filedialog.askopenfilename()

openfile = open(file, 'r')
cations = 0
anions = 0
NA = 0
CL = 0
for line in openfile:
    if 'ASP' in line or 'GLU' in line:
        anions = anions + 1
    if 'ARG' in line or 'LYS' in line:
        cations = cations + 1
    if 'NA+' in line:
        NA = NA + 1
    if 'CL-' in line:
        CL = CL + 1
print ('Cationic strength: ', cations)
print ('Anionic strength: ', anions)
print ('Net charge :', cations - anions)
print ('Na+ count: ' + str(NA))
print ('Cl- count: ' + str(CL))
      
openfile.close()
  