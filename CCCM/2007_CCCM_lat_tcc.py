# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will create a datasset of global total cloud cover with latitude.
Data is stored in a 2D array cccm_tcc_lat 
[:,0] = latitude
[:,1] = total cloud cover fraction
"""
import time
import numpy as np
import os
from pyhdf import SD
import matplotlib.pyplot as plt
import h5py

#---get latitude---#
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/') #Home PC
f = h5py.File('2007_CCCM_profile_variables.h5', 'r')

lat = f['lat'][:]

cf=[] # create a blank array to add cloud amount data

counter=0

# The directory where your HDF files are stored
os.chdir('E:/University/University/MSc/Models/Data/CCCM/2007')  # Home PC

#start = time.time()
# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the cloud cover data as a list. Since this is the 'cloud free' data, need to invert later
    cf = cf+f.select('Cloud area enhanced').get().tolist()

if len(lat) != len(cf):
    exit('Invalid sizes of lat and cf data')
    
cf = np.array(cf)

#Set the large 'fill values' in the data to nan before averaging        
cf[cf > 100] = np.nan

# Join the two lists as if they were two columns side by side, into a list of two elements each
combined = np.vstack((lat, cf)).T
#print ("combined")
#print (combined)

#print("get unique lats")
unique = np.unique(lat)
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print ("sorted")
#print (combined)

# Averages of (lat,cloud cover) empty array
averages_total = unique.size
cccm_tcc_lat = np.empty((averages_total,2),dtype=float)


# Current subtotal of current lat
subtotal = 0.0
# Current number of cloud cover entries in subtotal
number = 0
# Set the current lat to false
current_lat = None

# Iterate through all of the (lat,cloud cover) elements and subtotal the same lat values
i = 0
for item in combined:
    
    if np.isnan(item[1]):
        continue

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
        cccm_tcc_lat[i] = [current_lat, average]
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
#average = subtotal / number / 100
cccm_tcc_lat[i] = [current_lat, average]

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

plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tcc_lat[:,0],cccm_tcc_lat[:,1], 'b', label='CCCM')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

           
ax.set_ylabel('Cloud Fraction', color='r')
ax.set_xlabel('Latitude')

plt.title('Global Cloud Fraction vs Latitude - 2007')
plt.show()

os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM')
# specify path and file name to create 
with h5py.File('2007_CCCM_tcc_lat.h5', 'w') as p:
    p.create_dataset('tcc', data=cccm_tcc_lat)
    p.close()
