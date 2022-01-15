''# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:38:08 2017

Reads DSSP output from VMD and converts it into P(SS) vs. residue number.

@author: Erik Cook
"""

from tkinter import filedialog
import tkinter   
import numpy as np
import matplotlib.pyplot as plt

#correlation time for the electron nuclear interaction
tau = 0.0000000207
#Larmor frequency of the nuclear spin (proton, carbon, etc.) For proton, omega = 500,000,000Hz. For carbon, omega = 125,700,000Hz. 
H1omega = 498.186000
C13omega = 125.7
K = 1.23e-32        
t = 0.009
H1R2 = 0.0001
C13R2 = 0.00009
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
file = 'PRE Pdx1C C227S, S273C (CON).txt'   

data = np.loadtxt(file)

plt.figure(figsize=(10,6))
figure = plt.bar(data[:,0], data[:,1], color = "royalblue", edgecolor = "white", width = 1)


plt.tick_params(axis='both', which='major', labelsize=16)
plt.ylabel('$I_{para}/I_{dia}$', fontsize = 20)
plt.xlabel('Residue index', fontsize = 20)
plt.xlim(np.min(data[:,0]), np.max(data[:,0]))
plt.ylim(0,1.2)
C13R2SP = np.array([])
H1R2SP = np.array([])


aveDistance = np.loadtxt('N_012_AveDist_S273.txt')

H1R2SP = R2SP(tau, H1omega, K, aveDistance)
H1ratio = Iratio(H1R2, H1R2SP, t)
H1ratio = H1ratio - np.min(H1ratio)
H1ratio = H1ratio/np.max(H1ratio)

C13R2SP = R2SP(tau, C13omega, K, aveDistance)
C13ratio = Iratio(C13R2, C13R2SP, t)
C13ratio = C13ratio - np.min(C13ratio)
C13ratio = C13ratio/np.max(C13ratio)

residue = np.arange(203, 286, 1)

plt.plot(residue, C13ratio, color = 'black', linewidth = '2.5')    
plt.savefig('N_012_CON.pdf')
  