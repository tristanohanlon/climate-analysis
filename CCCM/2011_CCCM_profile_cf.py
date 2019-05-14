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
import h5py

#---get altitude in km---#
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/') #Home PC
f = h5py.File('2011_CCCM_profile_variables.h5', 'r')

lat = f['lat'][:]
alt_c = f['alt'][:]

#Empty lists
cf = []

# The directory where your HDF files are stored
#os.chdir('C:/Users/toha006/University/University/MSc/Models/Data/CCCM/2011')  # Uni Laptop
os.chdir('E:/University/University/MSc/Models/Data/CCCM/2011')  # Home PC

start = time.time()

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the cloud fraction data data which is a (lat, alt_c) array per file. Axis = 0 is added to for each file.  
    cf = cf+f.select('Cloud fraction profile').get().tolist() # 113 levels, #25535 values each referenced to latitude and longitude
    
end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()

cf = np.array(cf) # Convert cloud fraction list to a numpy array

#Set the large 'fill values' in the data to nan before averaging        
cf[cf > 100] = None

###############################################################################
#---Get Southern Ocean Data - Remember to get Latitude Data---#

# Join the two lists as if they were two columns side by side, into a list of two elements each
a = np.vstack(lat)
combined = np.hstack((a, cf))
#print ("combined")
#print (combined)
#Note: np.append(cf, lat, axis = 1) should work but will add the lat values to the end of the array

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print ("sorted")
#print (combined)

#Select latitudes over the southern ocean
combined = combined[combined[:,0]>=-70]
combined = combined[combined[:,0]<=-50]

#Split the combined array into just the cf data, eliminating the first coloumn of latitude
cf_so = combined[:,1:115] #only get values above 0 and above ground level
#alt_c = alt_c[0:112] #scale alt if necessary

###############################################################################

# Average the cloud fraction data over latitude and longitude for each altitude level        
cf = np.nanmean(cf, axis=0)
cf = np.vstack((alt_c, cf)).T

cf_so = np.nanmean(cf_so, axis=0)
cf_so = np.vstack((alt_c, cf_so)).T


#cf = cf[10:109] #only get values above 0 and above ground level
end = time.time()
print('Average data set creation took:', end - start, 's')

###########################################
#import h5py and save raw latitude data for iw and lw profile so use

os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM')
# specify path and file name to create 
with h5py.File('2011_CCCM_tcc_alt.h5', 'w') as p:
    p.create_dataset('cf', data=cf)
    p.create_dataset('cf_so', data=cf_so)
    p.close()