# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 17:55:00 2020

Calculates secondary structure of a trajectory using a PDB input file. Requires MDTraj.
Outputs a file containing DSSP values per residue for each model.

@author: St. Peppers
"""


import mdtraj as md
import time
import sys
import os

#USER DEFINED VARIABLES
#Construct file
extension = '_campari_trajs_traj.pdb'
runs = ['N_000', 'N_001', 'N_002', 'N_003', 'N_004', 'N_005', 'N_006', 'N_007', 'N_008', 'N_009', 'N_010', 'N_011', 'N_012', 'N_013', 'N_014', 'N_015']
#Add the extension to each run name.
filenames = []
for run in runs:
    filenames.append(run+extension)
    
    
#Start DSSP routine on each file sequentially.
for filename in filenames:
    file = open(filename, 'r')
    
    #Create output file to write DSSP calculations
    output = open(filename+'_DSSP.txt', 'w')
    output.close()
    
    #Create file to work from.
    workFileName = 'working.pdb'
    workFile = open(workFileName, 'w')
    
    #Model count
    modelsProcessed = 0
    
    # #Calculate the total number of models
    # totalModels = 0
    # for line in file:
    #     if "ENDMDL" in line.split() or "END" in line.split():
    #         totalModels += 1
    #     else:
    #         continue
    # sys.stdout.write("Performing DSSP calculation on " + str(totalModels) + " models.\n")
    
    sys.stdout.write('\nStarting DSSP calculations on ' + str(filename) + '...\n')
    
    for line in file:
        if "ENDMDL" in line.split() or "END" in line.split():
            workFile.close()
            #PROCESS workFile which contains one frame from the trajectory.
            traj = md.load_pdb(workFileName)
            dssp = md.compute_dssp(traj)
            output = open(filename+'_DSSP.txt', "a")
            output.write("\n")
            for i in dssp[0]:
                output.write(str(i) + " ")
            output.close()
            
            modelsProcessed += 1
            #FLUSH workFile
            workFile = open(workFileName, 'w')
            sys.stdout.write("\r" + "Models processed: " + str(modelsProcessed))
            sys.stdout.flush()
            
        elif "ATOM" in line.split():
            workFile.write(line)
        
        else:
            continue
    
    #Clean up workspace.
    sys.stdout.write('\nFinishing processing models. Cleaning workspace and processing output.\n')
    sys.stdout.write(str(modelsProcessed) + " model were processed using DSSP.")
    workFile.close()
    os.remove(workFileName)


    
