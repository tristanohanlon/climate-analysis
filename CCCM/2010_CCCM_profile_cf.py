# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will extract the cloud fraction profile against altitude. 
The cloud fraction content is averaged over all longitude and latitude at each altitude layer.
"""

import numpy as np
import os
from pyhdf import SD

#Empty lists
cf = []
alt_c = []

# The directory where your HDF files are stored
os.chdir('2010') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the altitude levels which corespond to the cloud fraction data. Only store the last file's values, already in km.
    alt_c = f.select('Layer center height profile (clouds and aerosol)').get().tolist() #113 levels for clouds and aerosols
    # Get the cloud fraction data data which is a (25535, 113) array per file. Axis = 0 is added to for each file.  
    cf = cf+f.select('Cloud fraction profile').get().tolist() # 113 levels, #25535 values each referenced to latitude and longitude

cf = np.array(cf) # Convert cloud fraction list to a numpy array

#Set the large 'fill values' in the data to 0 before averaging  
for index, item in np.ndenumerate(cf):
    if (item > 100): # 100% is the upper range limit of the cloud fraction data
        cf[index] = 0

# Average the cloud fraction data over latitude and longitude for each altitude level        
cf = np.mean(cf, axis=0)
