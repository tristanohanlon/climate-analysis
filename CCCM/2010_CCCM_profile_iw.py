# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will extract the global ice water content profile against altitude. 
The ice water content is averaged over all longitude and latitude at each altitude layer.
"""
import time
import numpy as np
import os
from pyhdf import SD

#Empty lists
alt = []
iw =[]
lat = []

start = time.time()

# The directory where your HDF files are stored
#os.chdir('C:/Users/toha006/University/University/MSc/Models/Data/CCCM/2010')  # Uni Laptop
os.chdir('E:/University/University/MSc/Models/Data/CCCM/2010')  # Home PC

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the altitude levels which corespond to the ice water content data. Only store the last file's values
    alt = f.select('Irradiance layer center height profile').get().tolist() #137 levels for ice and liquid water profiles
    # Get the ice water content data which is a (25535, 138) array per file. Axis = 0 is added to for each file.  
    iw = iw+f.select('Ice water content profile used').get().tolist() # 137 levels, 25536 values each referenced to latitude and longitude
    # Get the latitude data as a list
#    lat = lat+f.select('Colatitude of subsatellite point at surface at observation').get().tolist()
    
end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()

#lat[:] = [(round(v*2,0)/2-90)*-1 for v in lat]
#print("round lats")
#lat = np.array(lat)
alt = np.array(alt) / 1000 # Convert the altitude list to a numpy array and convery m to km
iw = np.array(iw) # Convert ice water content list to a numpy array

#Set the large 'fill values' in the data to nan before averaging
iw[iw > 20] = None        
#for index, item in np.ndenumerate(iw):
#    if (item > 20): # 20g/m^3 is the upper range limit of the ice water content data
#        iw[index] = 0     

# Join the two lists as if they were two columns side by side, into a list of two elements each
#combined = np.vstack((lat, iw)).T
#print ("combined")
#print (combined)

# Add a column for every additional column, -1 will sort by the first column
#combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print ("sorted")
#print (combined)

#Select latitudes over the southern ocean
#combined = combined[combined[:,0]>=-70]
#combined = combined[combined[:,0]<=-50]

#Split the combined array into just the iw data, eliminating the first coloumn of latitude
#iw = iw[:,1:138]

# Average the ice water content over latitude and longitude for each altitude level 
iw = np.nanmean(iw, axis=0)
iw = np.vstack((alt, iw)).T
end = time.time()
print('Average data set creation took:', end - start, 's')