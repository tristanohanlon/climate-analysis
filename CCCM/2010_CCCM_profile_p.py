# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 15:13:13 2019

@author: Tristan O'Hanlon

This will extract the global pressure profile against altitude. 
The pressure is averaged over all longitude and latitude at each altitude layer.
"""
import numpy as np
import os
from pyhdf import SD

#Empty lists
pressure = []
alt_t = []

# The directory where your HDF files are stored
os.chdir('2010') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the altitude levels which corespond to the pressure data. Only store the last file's values
    alt_t = f.select('Irradiance level height profile').get().tolist() #138 altitude levels pressure. 
    #Get the pressure data which is a (25535, 138) array per file. Axis = 0 is added to for each file.
    pressure = pressure+f.select('Pressure profile').get().tolist()

alt_t = np.array(alt_t) / 1000 # Convert the altitude list to a numpy array and convery m to km
pressure = np.array(pressure) # Convert pressure list to a numpy array

#Set the large 'fill values' in the data to 0 before averaging
for index, item in np.ndenumerate(pressure):
    if (item > 1100): # 1100hPa is the upper range limit of the pressure data
        pressure[index] = 0
# Average the pressure over latitude and longitude for each altitude level 
pressure = np.mean(pressure, axis=0)