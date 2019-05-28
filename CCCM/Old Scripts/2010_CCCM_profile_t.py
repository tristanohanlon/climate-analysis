# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 15:13:13 2019

@author: Tristan O'Hanlon

This will extract the global temperature profile against altitude. 
The temperature is averaged over all longitude and latitude at each altitude layer.
"""
import time
import numpy as np
import os
from pyhdf import SD

#Empty lists
temp=[]
#alt_t = []
#lat = []

start = time.time()

# The directory where your HDF files are stored
#os.chdir('C:/Users/toha006/University/University/MSc/Models/Data/CCCM/2010')  # Uni Laptop
os.chdir('E:/University/University/MSc/Models/Data/CCCM/2010')  # Home PC

# Load every file in the directory
for item in os.listdir(): 
    
    # Load the file
    f = SD.SD(item)
    # Get the altitude levels which corespond to the temperature data. Only store the last file's values
#    alt_t = f.select('Irradiance level height profile').get().tolist() #138 levels for temperature in (km)
    #Get the temperature data data which is a (25535, 138) array per file. Axis = 0 is added to for each file.
    temp = temp+f.select('Temperature profile').get().tolist()
    # Get the latitude data as a list
#    lat = lat+f.select('Colatitude of subsatellite point at surface at observation').get().tolist()
    
end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()

#lat[:] = [(round(v*2,0)/2-90)*-1 for v in lat]
#print("round lats")
#lat = np.array(lat)
#alt_t = np.array(alt_t)# Convert the altitude list to a numpy array
temp = np.array(temp) # Convert temperature list to a numpy array

#Set the large 'fill values' in the data to nan before averaging
temp[temp > 400] = None  
#for index, item in np.ndenumerate(temp):
#    if (item > 400): # 400K is the upper range limit of the temperature data
#        temp[index] = 0

# Join the two lists as if they were two columns side by side, into a list of two elements each
#lat = np.vstack(lat)
combined = np.hstack((lat, temp))
#print ("combined")
#print (combined)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print ("sorted")
#print (combined)

#Select latitudes over the southern ocean
combined = combined[combined[:,0]>=-70]
combined = combined[combined[:,0]<=-50]

#Split the combined array into just the temp data, eliminating the first coloumn of latitude
temp = temp[:,1:139]
#alt_t = alt_t[0:137] #scale alt if necessary

# Average the temperature over latitude and longitude for each altitude level     
temp = np.nanmean(temp, axis=0)
temp = np.vstack((alt_t, temp)).T

temp = temp[0:134] #to select alitude range
end = time.time()
print('Average data set creation took:', end - start, 's')