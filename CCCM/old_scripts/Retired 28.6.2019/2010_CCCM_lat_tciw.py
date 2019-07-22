# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will create a datasset of specific cloud ice water content (kgm^-2) with latitude.
Data is stored in a 2D array cccm_tciw_lat 
[:,0] = latitude
[:,1] = IWP
"""
import time
import numpy as np
from scipy import integrate
from sklearn.impute import SimpleImputer
import os
from pyhdf import SD
import matplotlib.pyplot as plt
import h5py

#---get altitude in km---#
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/') #Home PC
f = h5py.File('2010_CCCM_profile_variables.h5', 'r')

lat = f['lat'][:]
alt = f['alt'][:]

alt = alt*1000

#------------------------#

tciw = [] # create a blank array to add cloud ice water content data

os.chdir('E:/University/University/MSc/Models/Data/CCCM/2010')  # Home PC

start = time.time()
# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the cloud ice water content data as a list. (25536, 137) 'units': 'grams per cubic meter'
    tciw = tciw+f.select('Ice water content profile used').get().tolist() #same as profile plots
    
if len(lat) != len(tciw):
    exit('Invalid sizes of lat and tciw data')
    
end = time.time()
print('Importing data from files to lists took:', end - start, 's')

tciw = np.array(tciw)

#Set the large 'fill values' in the data to nan before averaging        
tciw[tciw > 20] = np.nan

####################

#fit all nan values to average

imp = SimpleImputer(missing_values=np.nan, strategy='mean')
imp.fit(np.transpose(tciw))  
a = imp.transform(np.transpose(tciw))
tciw = np.transpose(a)

####################
#computing the total cloud liquid water cloud content (LWP) kgm^-2

s_tciw = integrate.trapz(tciw, alt) # integrate across total altitude
a_tciw = -s_tciw / 1000 #convert g to kg

# Join the two lists as if they were two columns side by side, into a list of two elements each
combined = np.vstack((lat, a_tciw)).T
#print ("combined")
#print (combined)

#print("get unique lats")
unique = np.unique(lat)
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print ("sorted")
#print (combined)

# Averages of (lat, cloud ice water content) empty array
averages_total = unique.size
cccm_tciw_lat = np.empty((averages_total,2),dtype=float)

# Current subtotal of current lat
subtotal = 0.0
# Current number of cloud ice water content entries in subtotal
number = 0
# Set the current lat to false
current_lat = None

# Iterate through all of the (lat, cloud ice water content) elements and subtotal the same lat values
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
        average = subtotal / number
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
    subtotal+=item[1]
    
# Catch the last entry in the for loop
average = subtotal / number
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

 
plt.figure()
fig, ax = plt.subplots()
#ax.plot(cccm_tciw_frac_lat[:,0],cccm_tciw_frac_lat[:,1], 'g', label='CCCM')
ax.plot(cccm_tciw_lat[:,0],cccm_tciw_lat[:,1], 'g', label='CCCM')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

#ax.set_ylabel('Cloud Ice Water Content Fraction', color='r')           
ax.set_ylabel('Cloud Ice Water Content kg$m^{-2}$', color='r')
ax.set_xlabel('Latitude')

plt.title('IWP vs Latitude - 2010')
plt.show()

#import h5py

os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM')
# specify path and file name to create 
with h5py.File('2010_CCCM_tciw_lat.h5', 'w') as p:
    p.create_dataset('tciw', data=cccm_tciw_lat)
    p.close()
