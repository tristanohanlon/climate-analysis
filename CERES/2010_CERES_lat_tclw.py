# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will create a dataset of global liquid water fraction of total cloud cover with latitude.
Data is stored in a 2D array ceres_tclw_lat
[:,0] = latitude
[:,1] = liquid water fraction
"""

import numpy as np
import os
from pyhdf import SD

lw=[]
lat=[] # create a blank array to add latitude data
counter=0

# The directory where your HDF files are stored
os.chdir('2010/') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    
    # Get the latitude data as a list
    lat = lat+f.select('latitude').get().tolist()
    
    # Get the liquid water cloud amount data as a list. Zonal data is 5,180 so is averaged over longitude already.
    # The 5th index of axis = 0 gives the total liquid water cloud amount over latitude
    lw = lw+f.select('cld_amount_liq_zon').get()[-1,:].tolist()
    
if len(lat) != len(lw):
    exit('Invalid sizes of lat and cf data')
    
# Join the two lists as if they were two columns side by side, into a list of two elements each

lat = np.array(lat)
lw = np.array(lw)

combined_lw = np.vstack((lat, lw)).T

#print("unique lats")
unique = np.unique(lat)
#print(unique)

#print ("combined")
#print (combined)

#print ("sorted")
# Add a column for every additional column, -1 will sort by the first column
combined_lw = combined_lw[np.lexsort(np.transpose(combined_lw)[:-1])]
#print (combined)

# Averages of (lat,cloud cover) array
tclw = np.empty((unique.size,2),dtype=float)

# Current subtotal of current lat
subtotal = 0.0
# Current number of liquid water fraction entries in subtotal
number = 0
# Set the current lat to false
current_lat = None

# Iterate through all of the (lat,liquid water fraction) elements and subtotal the same lat values
i = 0

for item in combined_lw:

    if current_lat is None:
        current_lat = item[0];
    
    # If the lat is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_lat:
        
        # Find the average value
        average = subtotal / number / 100
        # Append the average
        tclw[i] = [current_lat, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current latitude
        current_lat = item[0]
        # Move to the next index in the averages array
        i+=1
    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number / 100
tclw[i] = [current_lat, average]
