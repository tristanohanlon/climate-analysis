# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will create a datasset of specific cloud liquid water content (kgm^-2) with latitude.
Data is stored in a 2D array cccm_tclw_lat 
[:,0] = latitude
[:,1] = LWP
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
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets') #Home PC
f = h5py.File('2006_CCCM_profile_variables.h5', 'r')

lat = f['lat'][:]
alt = f['alt'][:]
cff = f['cff'][:]

alt = alt*1000

#------------------------#

tclw = [] # create a blank array to add cloud liquid water content data

os.chdir('E:/University/University/MSc/Models/Data/CCCM/2006')  # Home PC

start = time.time()
# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the cloud liquid water content data as a list. (25536, 137) 'units': 'grams per cubic meter'
    tclw = tclw+f.select('Liquid water content profile used').get().tolist() #same as profile plots
    
if len(lat) != len(tclw):
    exit('Invalid sizes of lat and tclw data')
    
end = time.time()
print('Importing data from files to lists took:', end - start, 's')

tclw = np.array(tclw)

#Set the large 'fill values' in the data to nan before averaging        
tclw[tclw > 20] = np.nan

####################

#fit all nan values to average

imp = SimpleImputer(missing_values=np.nan, strategy='mean')
imp.fit(np.transpose(tclw))  
a = imp.transform(np.transpose(tclw))
tclw = np.transpose(a)

####################
#computing the total cloud liquid water cloud content (LWP) kgm^-2

s_tclw = integrate.trapz(tclw, alt) # integrate across total altitude
a_tclw = -s_tclw / 1000 #convert g to kg



# Join the two lists as if they were two columns side by side, into a list of two elements each
combined = np.vstack((lat, a_tclw)).T
#print ("combined")
#print (combined)

#print("get unique lats")
unique = np.unique(lat)
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print ("sorted")
#print (combined)

# Averages of (lat, cloud liquid water content) empty array
averages_total = unique.size
cccm_tclw_lat = np.empty((averages_total,2),dtype=float)

# Current subtotal of current lat
subtotal = 0.0
# Current number of cloud liquid water content entries in subtotal
number = 0
# Set the current lat to false
current_lat = None

# Iterate through all of the (lat, cloud liquid water content) elements and subtotal the same lat values
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
average = subtotal / number
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

cccm_tclw_lat = cccm_tclw_lat[:,1]*cff
cccm_tclw_lat = np.vstack((unique,cccm_tclw_lat)).T
 
plt.figure()
fig, ax = plt.subplots()
#ax.plot(cccm_tclw_frac_lat[:,0],cccm_tclw_frac_lat[:,1], 'g', label='CCCM')
ax.plot(cccm_tclw_lat[:,0],cccm_tclw_lat[:,1], 'g', label='CCCM')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

#ax.set_ylabel('Cloud liquid Water Content Fraction', color='r')           
ax.set_ylabel('Cloud liquid Water Content kg$m^{-2}$', color='r')
ax.set_xlabel('Latitude')

plt.title('LWP vs Latitude - 2006')
plt.show()

#import h5py

os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets')
# specify path and file name to create 
with h5py.File('2006_CCCM_tclw_lata.h5', 'w') as p:
    p.create_dataset('tclw', data=cccm_tclw_lat)
    p.close()
