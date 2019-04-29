# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 15:13:13 2019

@author: Tristan O'Hanlon

This will extract the global temperature profile against altitude. 
The temperature is averaged over all longitude and latitude at each altitude layer.
"""
import numpy as np
import os
from pyhdf import SD

#Empty lists
temp=[]
alt_t = []
counter=0

# The directory where your HDF files are stored
os.chdir('2010') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the altitude levels which corespond to the temperature data. Only store the last file's values
    alt_t = f.select('Irradiance level height profile').get().tolist() #138 levels for temperature
    #Get the temperature data data which is a (25535, 138) array per file. Axis = 0 is added to for each file.
    temp = temp+f.select('Temperature profile').get().tolist()


alt_t = np.array(alt_t) / 1000 # Convert the altitude list to a numpy array and convery m to km
temp = np.array(temp) # Convert temperature list to a numpy array

#Set the large 'fill values' in the data to 0 before averaging
for index, item in np.ndenumerate(temp):
    if (item > 400): # 400K is the upper range limit of the temperature data
        temp[index] = 0
        
# Average the temperature over latitude and longitude for each altitude level     
temp = np.mean(temp, axis=0)
