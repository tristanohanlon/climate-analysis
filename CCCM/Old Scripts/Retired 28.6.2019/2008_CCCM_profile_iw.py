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
import h5py

#Empty lists
iw =[]

#Import altitude and latitude data from reduced datasets

#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Uni Laptop
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/') #Home PC
p = h5py.File('2008_CCCM_profile_variables.h5', 'r')

alt = p['alt'][:]
lat = p['lat'][:]
air_density_g = p['air_density_g'][:]
air_density_so = p['air_density_so'][:]


###############################################################################

start = time.time()

# The directory where your HDF files are stored
os.chdir('E:/University/University/MSc/Models/Data/CCCM/2008')  # Home PC

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the ice water content (gm^-3) data which is a (lat, alt) array per file. Axis = 0 is added to for each file.  
    iw = iw+f.select('Ice water content profile used').get().tolist() # 137 levels, 25536 values each referenced to latitude and longitude
    
end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()

iw = np.array(iw) # Convert ice water content list to a numpy array

#Set the large 'fill values' in the data to nan before averaging
iw[iw > 20] = None        

###############################################################################
#---Get Southern Ocean Data - Remember to get Latitude Data---#

# Join the two lists as if they were two columns side by side, into a list of two elements each
a = np.vstack(lat)
combined = np.hstack((a, iw))
#print ("combined")
#print (combined)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print ("sorted")
#print (combined)

#Select latitudes over the southern ocean
combined = combined[combined[:,0]>=-70]
combined = combined[combined[:,0]<=-50]

#Split the combined array into just the iw data, eliminating the first coloumn of latitude
iw_so = combined[:,1:138]
#alt = alt[0:137] #scale alt if necessary


###############################################################################

# Average the ice water content over latitude and longitude for each altitude level 
iw = np.nanmean(iw, axis=0)
iw_so = np.nanmean(iw_so, axis=0)

###############################################################################

#Calculate density profile and divide IWC to get specific IWC at each profile level

iw = iw / air_density_g #kg/kg Global

iw_so = iw_so / air_density_so #kg/kg Southern Ocean


###############################################################################

iw = np.vstack((alt, iw)).T
iw_so = np.vstack((alt, iw_so)).T

end = time.time()
print('Average data set creation took:', end - start, 's')


###############################################################################
import h5py

os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM')
# specify path and file name to create 
with h5py.File('2008_CCCM_tciw_alt.h5', 'w') as p:
    p.create_dataset('iw', data=iw)
    p.create_dataset('iw_so', data=iw_so)
    p.close()