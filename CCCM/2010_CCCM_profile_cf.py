# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will extract the cloud fraction profile against altitude. 
The cloud fraction content is averaged over all longitude and latitude at each altitude layer.
"""
import time
import numpy as np
import os
from pyhdf import SD

#Empty lists
cf = []
alt_c = []
lat = []

# The directory where your HDF files are stored
#os.chdir('C:/Users/toha006/University/University/MSc/Models/Data/CCCM/2010')  # Uni Laptop
os.chdir('E:/University/University/MSc/Models/Data/CCCM/2010')  # Home PC

start = time.time()

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the altitude levels which corespond to the cloud fraction data. Only store the last file's values, already in km.
    alt_c = f.select('Layer center height profile (clouds and aerosol)').get().tolist() #113 levels for clouds and aerosols
    # Get the cloud fraction data data which is a (25535, 113) array per file. Axis = 0 is added to for each file.  
    cf = cf+f.select('Cloud fraction profile').get().tolist() # 113 levels, #25535 values each referenced to latitude and longitude
    # Get the latitude data as a list
    lat = lat+f.select('Colatitude of subsatellite point at surface at observation').get().tolist()
    
end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()

lat[:] = [(round(v*2,0)/2-90)*-1 for v in lat]
#print("round lats")
lat = np.array(lat)
cf = np.array(cf) # Convert cloud fraction list to a numpy array

#Set the large 'fill values' in the data to nan before averaging        
cf[cf > 100] = None
#for index, item in np.ndenumerate(cf):
#    if (item > 100): # 100% is the upper range limit of the cloud fraction data
#        cf[index] = 0

# Join the two lists as if they were two columns side by side, into a list of two elements each
combined = np.vstack((lat, cf)).T
#print ("combined")
#print (combined)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print ("sorted")
#print (combined)

#Select latitudes over the southern ocean
#combined = combined[combined[:,0]>=-70]
#combined = combined[combined[:,0]<=-50]

#Split the combined array into just the cf data, eliminating the first coloumn of latitude
#cf = cf[:,1:114]

# Average the cloud fraction data over latitude and longitude for each altitude level        
cf = np.nanmean(cf, axis=0)

end = time.time()
print('Average data set creation took:', end - start, 's')