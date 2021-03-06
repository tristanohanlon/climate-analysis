# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 15:13:13 2019

@author: Tristan O'Hanlon

This will extract the global pressure profile against altitude. 
The pressure is averaged over all longitude and latitude at each altitude layer.
"""
import time
import numpy as np
import os
from pyhdf import SD

#Empty lists
pressure = []
#alt_t = []
#lat = []

start = time.time()

# The directory where your HDF files are stored
#os.chdir('C:/Users/toha006/University/University/MSc/Models/Data/CCCM/2010')  # Uni Laptop
os.chdir('E:/University/University/MSc/Models/Data/CCCM/2010')  # Home PC

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the altitude levels which corespond to the pressure data. Only store the last file's values
#    alt_t = f.select('Irradiance level height profile').get().tolist() #138 altitude levels pressure. 
    #Get the pressure data which is a (25535, 138) array per file. Axis = 0 is added to for each file.
    pressure = pressure+f.select('Pressure profile').get().tolist()
    # Get the latitude data as a list
#    lat = lat+f.select('Colatitude of subsatellite point at surface at observation').get().tolist()
    
end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()

#lat[:] = [(round(v*2,0)/2-90)*-1 for v in lat]
#print("round lats")
#lat = np.array(lat)
#alt_t = np.array(alt_t) / 1000 # Convert the altitude list to a numpy array and convery m to km
pressure = np.array(pressure) # Convert pressure list to a numpy array

#Set the large 'fill values' in the data to nan before averaging
pressure[pressure > 1100] = None   
#for index, item in np.ndenumerate(pressure):
#    if (item > 1100): # 1100hPa is the upper range limit of the pressure data
#        pressure[index] = 0

# Join the two lists as if they were two columns side by side, into a list of two elements each
combined = np.hstack((lat, pressure))
#print ("combined")
#print (combined)

# Add a column for every additional column, -1 will sort by the first column
#combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print ("sorted")
#print (combined)

#Select latitudes over the southern ocean
combined = combined[combined[:,0]>=-70]
combined = combined[combined[:,0]<=-50]

#Split the combined array into just the iw data, eliminating the first coloumn of latitude
pressure = pressure[:,1:139]
#alt_t = alt_t[0:137] #scale alt if necessary

# Average the pressure over latitude and longitude for each altitude level 
pressure = np.nanmean(pressure, axis=0)
pressure = np.vstack((alt_t, pressure)).T
pressure = pressure[0:134]
end = time.time()
print('Average data set creation took:', end - start, 's')