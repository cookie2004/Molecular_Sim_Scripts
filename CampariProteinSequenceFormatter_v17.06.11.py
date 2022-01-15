# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 18:01:12 2017

Input a FASTA formatted single-letter sequence for a protein, and the Campari-formatted SEQ file is exported.

@author: St. Peppers
"""

seq = input('Input the single-letter sequence of your protein. \n')

output = []

for aa in seq:
    if aa.lower() == 'a':
        output.append('ALA')
    if aa.lower() == 'c':
        output.append('CYS')
    if aa.lower() == 'd':
        output.append('ASP')
    if aa.lower() == 'e':
        output.append('GLU')
    if aa.lower() == 'f':
        output.append('PHE')
    if aa.lower() == 'g':
        output.append('GLY')
    if aa.lower() == 'h':
        his = input('HIE or HID?\n')
        if his.lower() == 'hie':
            output.append('HIE')
        if his.lower() == 'hid':
            output.append('HID')
    if aa.lower() == 'i':
        output.append('ILE')
    if aa.lower() == 'k':
        output.append('LYS')
    if aa.lower() == 'l':
        output.append('LEU')
    if aa.lower() == 'm':
        output.append('MET')
    if aa.lower() == 'n':
        output.append('ASN')
    if aa.lower() == 'p':
        output.append('PRO')
    if aa.lower() == 'q':
        output.append('GLN')
    if aa.lower() == 'r':
        output.append('ARG')
    if aa.lower() == 's':
        output.append('SER')
    if aa.lower() == 't':
        output.append('THR')
    if aa.lower() == 'v':
        output.append('VAL')
    if aa.lower() == 'w':
        output.append('TRP')
    if aa.lower() == 'y':
        output.append('TYR')
       
outputfile = open("outputsequence.txt", 'w')  

for i in output:
    outputfile.write(str(i))
    outputfile.write("\n")     
    
outputfile.close()