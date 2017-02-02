# -*- coding: utf-8 -*-
"""
Plot the data amassed from HyPerCarlo
"""
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re

"""
Change font size, was way too small
"""
params = {'font.size':20,
          'legend.fontsize':16}
          
plt.rcParams.update(params)

path = 'data/DipoleDipole/'
runData = dict()
for infile in glob.glob( os.path.join(path, 'run*/Output.dat') ):
    print "current file: " + infile
    
    for line in open(infile, 'r'):
        if re.search('Nx', line):
            runInfo=re.sub(r'[\W_]+', '', line)
    runData[runInfo] = np.loadtxt(infile)
    
print runData

plt.figure(1)
plt.xlabel('T')
plt.ylabel('E/N (eV)')
for run in runData:
    plt.plot(runData[run][:,0], runData[run][:,1], label=run, marker='^', linestyle='--')
plt.legend(loc='lower right', numpoints=1)
plt.tight_layout()
plt.savefig('E_vs_T.pdf')
    
plt.figure(2)

plt.xlabel('T')
plt.ylabel('M/N')
for run in runData:
    plt.plot(runData[run][:,0], runData[run][:,3], label=run, marker='^', linestyle='--')
plt.legend(loc='upper right', numpoints=1)
plt.tight_layout()
plt.savefig('P_vs_T.pdf')

cvmaxT=0.0
cvmax=0.0
plt.figure(3)
plt.xlabel('T')
plt.ylabel(r'C$_v$/N')
for run in runData:
    if max(runData[run][:,5])>cvmax:
        cvmax=max(runData[run][:,5])
        cvmaxT=runData[run][(runData[run][:,5].argmax()),0]
    plt.plot(runData[run][:,0], runData[run][:,5], label=run, marker='^', linestyle='--')
plt.axvline(x=cvmaxT,ymin=0,ymax=1,linestyle='--', color='red', label=r'$C_v^{max}$ (T='+str(int(cvmaxT)) +')')
plt.legend(loc='upper right', numpoints=1)
plt.tight_layout()
plt.savefig('Cv_vs_T.pdf')

chimaxT=0.0
chimax=0.0

plt.figure(4)
plt.xlabel('T (K)')
plt.ylabel(r'$\chi$/N')
for run in runData:
    if max(runData[run][:,6])>chimax:
        chimax=max(runData[run][:,6])
        chimaxT=runData[run][(runData[run][:,6].argmax()),0]
    plt.plot(runData[run][:,0], runData[run][:,6], label=run, marker='^', linestyle='--')
plt.axvline(x=chimaxT,ymin=0,ymax=1,linestyle='--', color='red', label=r'$\chi^{max}$ (T='+str(int(chimaxT)) +')')
plt.legend(loc='upper right', numpoints=1)
plt.tight_layout()
plt.savefig('Chi_vs_T.pdf')
