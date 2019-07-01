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
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets') #Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets') #Uni laptop
#os.chdir('//synthesis/E/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets')

f = h5py.File('2006_CCCM_profile_variables.h5', 'r')

#lat = f['lat'][:]
alt = f['alt'][:]
air_density_g = f['air_density_g'][:]

#------------------------#

tciw = [] # create a blank array to add cloud ice water content data
lat = []
os.chdir('E:/University/University/MSc/Models/Data/CCCM/')  # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/Data/CCCM/Test') #Uni laptop
#os.chdir('//synthesis/E/University/University/MSc/Models/Data/CCCM/2006') #Uni laptop

start = time.time()
# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the cloud ice water content data as a list. (lat, alt) 'units': 'grams per cubic meter'
    tciw = tciw+f.select('Ice water content profile used').get().tolist()
    lat = lat+f.select('Colatitude of subsatellite point at surface at observation').get().tolist()
    
if len(lat) != len(tciw):
    exit('Invalid sizes of lat and tciw data')
    
end = time.time()
print('Importing data from files to lists took:', end - start, 's')


start = time.time()

tciw = np.array(tciw)
      
tciw[tciw > 20] = np.nan #Set the large 'fill values' in the data to nan before averaging  
lat[:] = [(round(v*2,0)/2-90)*-1 for v in lat] #convert colatitude to latitude and raound to nearest 0.5 degrees

####################

#fit all nan values to average
"""
imp = SimpleImputer(missing_values=np.nan, strategy='mean')
imp.fit(np.transpose(tciw))  
a = imp.transform(np.transpose(tciw))
tciw = np.transpose(a)
"""
####################

#a_tciw = tciw / 1000 #convert g to kg
a_lat = np.vstack(lat)

# Join the two lists as if they were two columns side by side, into a list of two elements each
combined = np.hstack((a_lat, tciw))
#print ("combined")
#print (combined)

#print("get unique lats")
unique = np.unique(lat)
#print(unique)

# Add a column for every additional column, -1 will sort by the first column

combined.sort(axis=0) #sort the combined array by latitude
#print ("sorted")
#print (combined)

#import h5py
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets')
#os.chdir('//synthesis/E/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets')
# specify path and file name to create 
with h5py.File('lw_combined.h5', 'w') as p:
    p.create_dataset('iw_alt_lat', data=cccm_tciw_lat)
    p.close()

# total bin size of latitude entires empty list
count_bin = np.zeros([329])
cccm_tciw_lat = np.empty((1,138),dtype=float)

# Current subtotal of current lat
subtotal = np.zeros([1, 138])
# Set the current lat to false
current_lat = None

end = time.time()
print('Getting loop data ready:', end - start, 's')


count = 0
i = 0
start = time.time()

# Iterate through all of the (lat, cloud ice water content) elements and subtotal the same lat values
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
        average = np.array([np.nansum(subtotal[1:], axis = 0)]) # exclude first row with zeros
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
        cccm_tciw_lat = np.concatenate((cccm_tciw_lat, average)) #stack new row
        # Reset the subtotal
        subtotal = np.zeros([1, 138])
        # Set the current latitude
        current_lat = item[0]
        # Move to the next index in the averages array
        count_bin[i] = count
        count = 0
        i += 1
    # Add the next value to the subtotal
    item = np.array([item]) # make the item list an array (1,138)
    subtotal=np.concatenate((subtotal, item)) #stack the item on the bottom of the subtotal array
    count += 1
   
    
end = time.time()
print('Computing average data:', end - start, 's')    
# Catch the last entry in the for loop
average = np.array([np.nansum(subtotal[1:], axis = 0)])
cccm_tciw_lat = np.concatenate((cccm_tciw_lat, average))
count_bin[i] = count

cccm_tciw_lat = cccm_tciw_lat[1:]

cccm_tciw_lat = cccm_tciw_lat / np.vstack(count_bin)
cccm_tciw_lat = cccm_tciw_lat[:,1:]

cccm_tciw_lat = cccm_tciw_lat / air_density_g


"""
print ("averages")
# Iterate through all of the (lat,cloud ice water content) elements
for item in averages:
    print("[", end='')
    print(item[0], end='')
    print(", ", end='')
    print(item[1], end='')
    print("]\n", end='')


cccm_tciw_lat = cccm_tciw_lat[:,1]*cff
cccm_tciw_lat = np.vstack((unique,cccm_tciw_lat)).T

plt.figure()
fig, ax = plt.subplots()
#ax.plot(cccm_tciw_frac_lat[:,0],cccm_tciw_frac_lat[:,1], 'g', label='CCCM')
ax.contourf(unique, alt[35:], np.transpose(cccm_tciw_lat[35:]))

#ax.set_ylabel('Cloud Ice Water Content Fraction', color='r')           
ax.set_ylabel('Altitude (km)', color='r')
ax.set_xlabel('Latitude')

plt.title('Ice Water Content - 2006')
plt.show()
"""

#import h5py
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets')
#os.chdir('//synthesis/E/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets')
# specify path and file name to create 
with h5py.File('CCCM_tciw_lat_alt.h5', 'w') as p:
    p.create_dataset('iw_alt_lat', data=cccm_tciw_lat)
    p.close()
