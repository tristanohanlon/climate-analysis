# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will create a dataset of global total cloud cover with latitude.
Data is stored in a 2D array ceres_tcc_lat
[:,0] = latitude
[:,1] = total cloud cover fraction
"""

import numpy as np
import os
from pyhdf import SD

cf=[] # create a blank array to add cloud amount data
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
    
    # Get the cloud cover data as a list. Zonal data is 5,180 so is averaged over longitude already. 
    # The 5th index of axis = 0 gives the total cloud cover over latitude
    cf = cf+f.select('cld_amount_zon').get()[-1,:].tolist()


if len(lat) != len(cf):
    exit('Invalid sizes of lat and cf data')


# Join the two lists as if they were two columns side by side, into a list of two elements each
lat = np.array(lat)
cf = np.array(cf)
combined = np.vstack((lat, cf)).T

#print("unique lats")
unique = np.unique(lat)
#print(unique)

#print ("combined")
#print (combined)

#print ("sorted")
# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print (combined)

# Averages of (lat,cloud cover) array
ceres_tcc_lat = np.empty((unique.size,2),dtype=float)

# Current subtotal of current lat
subtotal = 0.0
# Current number of cloud cover entries in subtotal
number = 0
# Set the current lat to false
current_lat = None

# Iterate through all of the (lat,cloud cover) elements and subtotal the same lat values
i = 0
for item in combined:

    if current_lat is None:
        """
        print("setting current_lat to item[0]")
        print("(current_lat == item[0]) = ", end='')
        print(current_lat == item[0]) 
        """
        current_lat = item[0];
    
    # If the lat is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_lat:
        
        # Find the average value
        average = subtotal / number / 100
        """
        print("--------")
        print("lat: ", end='')
        print(current_lat, end='')
        print(", avg: ", end='')
        print(average, end='')
        print(", subtotal: ", end='')
        print(subtotal, end='')
        print(", number: ", end='')
        print(number)
        """
        # Append the average
        ceres_tcc_lat[i] = [current_lat, average]
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
ceres_tcc_lat[i] = [current_lat, average]


"""
print ("averages")
# Iterate through all of the (lat,cloud cover) elements
for item in averages:
    print("[", end='')
    print(item[0], end='')
    print(", ", end='')
    print(item[1], end='')
    print("]\n", end='')
"""
