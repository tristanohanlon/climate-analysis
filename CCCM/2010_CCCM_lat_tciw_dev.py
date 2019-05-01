# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will create a datasset of global total cloud ice water content with latitude.
Data is stored in a 2D array cccm_tciw_lat 
[:,0] = latitude
[:,1] = total cloud ice water content
"""

import numpy as np
import os
from pyhdf import SD
import matplotlib.pyplot as plt

lat = []
tciw = []

#os.chdir('C:/Users/toha006/University/University/MSc/Models/Data/CCCM/2010')  # Uni Laptop
os.chdir('E:/University/University/MSc/Models/Data/CCCM/Test') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    
    
    # Get the latitude data as a list
    lat = lat+f.select('Colatitude of subsatellite point at surface at observation').get().tolist()
    
    # Get the cloud ice water content data as a list. (25536, 113) 'units': 'grams per cubic meter'
    tciw = tciw+f.select('Mean CloudSat radar only ice water content').get().tolist()

if len(lat) != len(tciw):
    exit('Invalid sizes of lat and tciw data')

# Join the two lists as if they were two columns side by side, into a list of two elements each
#print("round lats")
lat[:] = [(round(v*2,0)/2-90)*-1 for v in lat]

lat = np.array(lat)
tciw = np.array(tciw)

#Set the large 'fill values' in the data to 0 before averaging        
tciw[tciw > 10.0] = None
tciw = np.nanmean(tciw, axis=1)
combined = np.vstack((lat, tciw)).T

#print("get unique lats")
unique = np.unique(lat)
#print(unique)

#print ("combined")
#print (combined)

#print ("sorted")
# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print (combined)

# Averages of (lat,cloud ice water content) array
averages_total = unique.size
cccm_tciw_lat = np.empty((averages_total,2),dtype=float)

i = 0
     

# Current subtotal of current lat
subtotal = 0.0
# Current number of cloud ice water content entries in subtotal
number = 0
# Set the current lat to false
current_lat = None

# Iterate through all of the (lat,cloud ice water content) elements and subtotal the same lat values

for item in combined:

    if current_lat is None:
        """
        print("setting current_lat to item[0]")
        print("(current_lat == item[0]) = ", end='')
        print(current_lat == item[0]) 
        """
        current_lat = item[0];
    
    # If the lat is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_lat and item[1] != None:
        
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
        cccm_tciw_lat[i] = [current_lat, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current latitude
        current_lat = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal=np.nansum(item[1])
    
# Catch the last entry in the for loop
average = subtotal / number / 100
cccm_tciw_lat[i] = [current_lat, average]


"""
print ("averages")
# Iterate through all of the (lat,cloud ice water content) elements
for item in averages:
    print("[", end='')
    print(item[0], end='')
    print(", ", end='')
    print(item[1], end='')
    print("]\n", end='')
"""

#cccm_tciw_lat[cccm_tciw_lat==0] = None
   
plt.figure()
plt.plot(cccm_tciw_lat[:,0],cccm_tciw_lat[:,1])
plt.show()
