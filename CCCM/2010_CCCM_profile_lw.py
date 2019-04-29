# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will extract the global liquid water content profile against altitude. 
The liquid water content is averaged over all longitude and latitude at each altitude layer.
"""

import numpy as np
import os
from pyhdf import SD

#Empty lists
lw = []
alt = []

# The directory where your HDF files are stored
os.chdir('2010') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the altitude levels which corespond to the liquid water content data. Only store the last file's values
    alt = f.select('Irradiance layer center height profile').get().tolist() #137 levels for ice and liquid water profiles
    # Get the liquid water content data data which is a (25535, 138) array per file. Axis = 0 is added to for each file.  
    lw = lw+f.select('Liquid water content profile used').get().tolist() # 137 levels, 25536 values each referenced to latitude and longitude

alt = np.array(alt) / 1000 # Convert the altitude list to a numpy array and convery m to km
lw = np.array(lw) # Convert liquid water content list to a numpy array

#Set the large 'fill values' in the data to 0 before averaging        
for index, item in np.ndenumerate(lw):
    if (item > 20): # 20g/m^3 is the upper range limit of the liquid water content data
        lw[index] = 0
   
# Average the liquid water content over latitude and longitude for each altitude level 
lw = np.mean(lw, axis=0)