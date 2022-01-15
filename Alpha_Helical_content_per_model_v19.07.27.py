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
content = []
models = 0
file = "HISPEP19_1MGdn_Gdn1_2_SS.txt"
output = '1MGdn.txt'
numberOfRes = 0
data = open(file,'r')
stuff = np.loadtxt(file, skiprows=9, dtype=str)

for line in data:
    line = line.split(' ')
    if "NUM_ITEMS=" in line:
        numberOfRes = int(line[2])
    else:
        None
tempCont = []
for i in stuff:
    tempCont.append(i)
    
    if int(i[0]) == numberOfRes:
        #print (tempCont)
        a = 0
        for y in tempCont:
            if y[4] == 'H':
                a = a + 1
            else:
                None
        content.append(a/numberOfRes)
        #print (a/numberOfRes)
        tempCont = []
        models = models + 1
        print ('Models procssed: ' + str(models))
    else:
        None

plt.hist(content, 50)

#Export data.
out = open(output, 'w')
for i in content:
    out.write(str(i) + '\n')
out.close()