# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will extract the global liquid water content profile against altitude. 
The liquid water content is averaged over all longitude and latitude at each altitude layer.
"""
import time
import numpy as np
import os
from pyhdf import SD
import h5py

#Empty lists
lw =[]

#Import altitude and latitude data from reduced datasets

os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/') #Home PC
p = h5py.File('2007_CCCM_profile_variables.h5', 'r')

alt = p['alt'][:]
lat = p['lat'][:]
air_density_g = p['air_density_g'][:]
air_density_so = p['air_density_so'][:]


###############################################################################

start = time.time()

# The directory where your HDF files are stored
os.chdir('E:/University/University/MSc/Models/Data/CCCM/2007')  # Home PC

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the liquid water content (gm^-3) data which is a (lat, alt) array per file. Axis = 0 is added to for each file.  
    lw = lw+f.select('Liquid water content profile used').get().tolist() # 137 levels, 25536 values each referenced to latitude and longitude
    
end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()

lw = np.array(lw) # Convert liquid water content list to a numpy array

#Set the large 'fill values' in the data to nan before averaging
lw[lw > 20] = None        

###############################################################################
#---Get Southern Ocean Data - Remember to get Latitude Data---#

# Join the two lists as if they were two columns side by side, into a list of two elements each
a = np.vstack(lat)
combined = np.hstack((a, lw))
#print ("combined")
#print (combined)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print ("sorted")
#print (combined)

#Select latitudes over the southern ocean
combined = combined[combined[:,0]>=-70]
combined = combined[combined[:,0]<=-50]

#Split the combined array into just the lw data, eliminating the first coloumn of latitude
lw_so = combined[:,1:138]
#alt = alt[0:137] #scale alt if necessary


###############################################################################

# Average the liquid water content over latitude and longitude for each altitude level 
lw = np.nanmean(lw, axis=0)
lw_so = np.nanmean(lw_so, axis=0)

###############################################################################

#Calculate density profile and divide LWC to get specific lwC at each profile level

lw = lw / air_density_g #kg/kg Global

lw_so = lw_so / air_density_so #kg/kg Southern Ocean


###############################################################################

lw = np.vstack((alt, lw)).T
lw_so = np.vstack((alt, lw_so)).T

end = time.time()
print('Average data set creation took:', end - start, 's')


###############################################################################
import h5py

os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM')
# specify path and file name to create 
with h5py.File('2007_CCCM_tclw_alt.h5', 'w') as p:
    p.create_dataset('lw', data=lw)
    p.create_dataset('lw_so', data=lw_so)
    p.close()