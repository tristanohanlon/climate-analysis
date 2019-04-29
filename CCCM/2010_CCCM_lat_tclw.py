# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will create a datasset of global total cloud liquid water content with latitude.
Data is stored in a 2D array cccm_tclw_lat 
[:,0] = latitude
[:,1] = total cloud liquid water content
"""

import numpy as np
import os
from pyhdf import SD
import matplotlib.pyplot as plt

tclw=[] # create a blank array to add cloud liquid water content data
lat=[] # create a blank array to add latitude data
counter=0

# The directory where your HDF files are stored
os.chdir('2010/') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    
    
    # Get the latitude data as a list
    lat = lat+f.select('Colatitude of subsatellite point at surface at observation').get().tolist()
    
    # Get the cloud liquid water content data as a list. (25536, 113) 'units': 'grams per cubic meter'
    tclw = tclw+f.select('Mean CloudSat radar only liquid water content').get().tolist()

if len(lat) != len(cf):
    exit('Invalid sizes of lat and cf data')

# Join the two lists as if they were two columns side by side, into a list of two elements each
#print("round lats")
lat[:] = [(round(v*2,0)/2-90)*-1 for v in lat]

lat = np.array(lat)
tclw = np.array(tclw)

#Set the large 'fill values' in the data to 0 before averaging        
for index, item in np.ndenumerate(tclw):
    if (item >15): # 20g/m^3 is the upper range limit of the ice water content data
        tclw[index] = 0     

# Average the ice water content over latitude and longitude for each altitude level 
tclw = np.mean(tclw, axis=1)

combined = np.vstack((lat, tclw)).T

#print("get unique lats")
unique = np.unique(lat)
#print(unique)

#print ("combined")
#print (combined)

#print ("sorted")
# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print (combined)

# Averages of (lat,cloud liquid water content) array
averages_total = unique.size
cccm_tclw_lat = np.empty((averages_total,2),dtype=float)

# Current subtotal of current lat
subtotal = 0.0
# Current number of cloud liquid water content entries in subtotal
number = 0
# Set the current lat to false
current_lat = None

# Iterate through all of the (lat,cloud liquid water content) elements and subtotal the same lat values
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
        
        # Find the average value.
        average = subtotal / number / 100
        #average = subtotal / number / 100
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
        cccm_tclw_lat[i] = [current_lat, average]
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
cccm_tclw_lat[i] = [current_lat, average]


"""
print ("averages")
# Iterate through all of the (lat,cloud liquid water content) elements
for item in averages:
    print("[", end='')
    print(item[0], end='')
    print(", ", end='')
    print(item[1], end='')
    print("]\n", end='')
"""
   
plt.figure()
plt.plot(cccm_tclw_lat[:,0],cccm_tclw_lat[:,1])
plt.show()