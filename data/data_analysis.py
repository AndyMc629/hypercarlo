# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import glob
import os
import matplotlib.pyplot as plt

if __name__ == "__main__":
    
    path = r'/home/andrew/dev/GitHub/hypercarlo/data/run_24_01_2017_2119'                     # use your path
    
    runFiles = glob.glob(os.path.join(path, "Run*.dat"))     # advisable to use os.path.join as this makes concatenation OS independent
    equilFiles = glob.glob(os.path.join(path, "Equilibration*.dat"))
    
    dfFromRunFiles = (pd.read_csv(f, skiprows=2, sep=' ') for f in runFiles)
    dfFromEquilFiles = (pd.read_csv(f, skiprows=1, sep=' ') for f in equilFiles)
 
    
    dfRun = pd.concat(dfFromRunFiles, ignore_index=True)
    dfEquil = pd.concat(dfFromEquilFiles, ignore_index=True) 
    # doesn't create a list, nor does it append to one
    
    