# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will extract the cloud fraction profile against altitude. 
The global cloud fraction is already averaged over all longitude and latitude at each altitude layer.
[:,0] = altitude
[:,1] = cloud fraction
"""

import numpy as np
import os
from pyhdf import SD

cf=[] # create a blank array to add cloud amount data
alt=[] # create a blank array to add altitude data
counter=0

# The directory where your HDF files are stored
os.chdir('2010/') 

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    
    # Get the altitude data as a list
    alt = alt+f.select('cld_eff_hgt_glob').get()[:-1].tolist()
    
    # Get the cloud fraction data as a list, excluding the last value which is the total for all altitudes
    cf = cf+f.select('cld_amount_glob').get()[:-1].tolist()

if len(alt) != len(cf):
    exit('Invalid sizes of alt and cf data')


# Join the two lists as if they were two columns side by side, into a list of two elements each
alt = np.array(alt)
cf = np.array(cf)
combined = np.vstack((alt, cf)).T

#print("unique alts")
unique = np.unique(alt)
#print(unique)

#print ("combined")
#print (combined)

#print ("sorted")
# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print (combined)

# Averages of (alt,cloud fraction) array
ceres_tcc_alt = np.empty((unique.size,2),dtype=float)

# Current subtotal of current alt
subtotal = 0.0
# Current number of cloud fraction entries in subtotal
number = 0
# Set the current alt to false
current_alt = None

# Iterate through all of the (alt,cloud fraction) elements and subtotal the same alt values
i = 0
for item in combined:

    if current_alt is None:
        """
        print("setting current_alt to item[0]")
        print("(current_alt == item[0]) = ", end='')
        print(current_alt == item[0]) 
        """
        current_alt = item[0];
    
    # If the alt is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_alt:
        
        # Find the average value
        average = subtotal / number / 100
        """
        print("--------")
        print("lat: ", end='')
        print(current_alt, end='')
        print(", avg: ", end='')
        print(average, end='')
        print(", subtotal: ", end='')
        print(subtotal, end='')
        print(", number: ", end='')
        print(number)
        """
        # Append the average
        ceres_tcc_alt[i] = [current_alt, average]
        # Reset the subtotal
        subtotal = 0.0
        number = 0
        # Set the current altitude
        current_alt = item[0]
        # Move to the next index in the averages array
        i+=1

    # Add the next value to the subtotal
    number+=1
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number / 100
ceres_tcc_alt[i] = [current_alt, average]


"""
print ("averages")
# Iterate through all of the (alt,cloud fraction) elements
for item in averages:
    print("[", end='')
    print(item[0], end='')
    print(", ", end='')
    print(item[1], end='')
    print("]\n", end='')
"""
   


#plt.figure()
#plt.plot(ceres_tcc_lat[:,0],ceres_tcc_lat[:,1])
#plt.show()