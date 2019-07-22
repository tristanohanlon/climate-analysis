# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will create a datasset of global total cloud ice water content with latitude.
Data is stored in a 2D array cccm_tciw_lat 
[:,0] = latitude
[:,1] = total cloud ice water content
"""
import time
import numpy as np
import os
from pyhdf import SD
import matplotlib.pyplot as plt

lat = [] # create a blank array to add latitude data
tciw = [] # create a blank array to add cloud ice water content data

#os.chdir('C:/Users/toha006/University/University/MSc/Models/Data/CCCM/2010')  # Uni Laptop
os.chdir('E:/University/University/MSc/Models/Data/CCCM/2010')  # Home PC
#os.chdir('../Data/CCCM')  # Laptop

start = time.time()
# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
     
    # Get the latitude data as a list
    lat = lat+f.select('Colatitude of subsatellite point at surface at observation').get().tolist()
    # Get the cloud ice water content data as a list. (25536, 137) 'units': 'grams per cubic meter'
    tciw = tciw+f.select('Ice water content profile used').get().tolist() #same as profile plots
    
if len(lat) != len(tciw):
    exit('Invalid sizes of lat and tciw data')
    
end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()

lat[:] = [(round(v*2,0)/2-90)*-1 for v in lat]
#print("round lats")

lat = np.array(lat)
tciw = np.array(tciw)

#Set the large 'fill values' in the data to nan before averaging        
tciw[tciw > 20] = np.nan

tciw = tciw / air_density

# Average the ice water content over latitude and longitude for each altitude level 
tciw = np.nansum(tciw, axis=1)

# Join the two lists as if they were two columns side by side, into a list of two elements each
combined = np.vstack((lat, tciw)).T
#print ("combined")
#print (combined)

#print("get unique lats")
unique = np.unique(lat)
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print ("sorted")
#print (combined)

# Averages of (lat,cloud ice water content) empty array
averages_total = unique.size
cccm_tciw_lat = np.empty((averages_total,2),dtype=float)

#end = time.time()
#print('Create arrays and combined array took:', end - start, 's')

start = time.time()

# Current subtotal of current lat
subtotal = 0.0
# Current number of cloud ice water content entries in subtotal
number = 0
# Set the current lat to false
current_lat = None

# Iterate through all of the (lat,cloud ice water content) elements and subtotal the same lat values
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

end = time.time()
print('Average data set creation took:', end - start, 's')
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

# Make the data a fraction of the total ice and liquid contribution - multiplied by the total cloud fraction
# cccm_tciw_frac_lat =  cccm_tciw_lat
# cccm_tciw_frac_lat[:,1] = (cccm_tciw_lat[:,1]  / (cccm_tciw_lat[:,1] + cccm_tclw_lat[:,1])) * cccm_tciw_lat[:,1] + cccm_tcc_lat[:,1]

#Select latitudes over the southern ocean
#cccm_tciw_lat = cccm_tciw_lat[cccm_tciw_lat[:,0]>=-70]
#cccm_tciw_lat = cccm_tciw_lat[cccm_tciw_lat[:,0]<=-50]
   
plt.figure()
fig, ax = plt.subplots()
#ax.plot(cccm_tciw_frac_lat[:,0],cccm_tciw_frac_lat[:,1], 'g', label='CCCM')
ax.plot(cccm_tciw_lat[:,0],cccm_tciw_lat[:,1], 'g', label='CCCM')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

#ax.set_ylabel('Cloud Ice Water Content Fraction', color='r')           
ax.set_ylabel('Specific Cloud Ice Water Content kg/kg', color='r')
ax.set_xlabel('Latitude')

plt.title('Cloud Ice Water Content vs Latitude - 2010')
plt.show()
"""
import h5py

os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM')
# specify path and file name to create 
with h5py.File('2010_CCCM_SO_profile_data_test1.h5', 'w') as p:
    pcf = p.create_dataset('Cloud Fraction', data=cccm_tcc_lat)
    pcf.attrs['title'] = "Cloud Fraction"
    plw = p.create_dataset('Liquid Water Content', data=cccm_tclw_lat)
    plw.attrs['title'] = "Cloud Liquid Water Content"
    plw.attrs['units'] = "gm^-3" 
    piw = p.create_dataset('Ice Water Content', data=cccm_tciw_lat)
    piw.attrs['title'] = "Cloud Ice Water Content"
    piw.attrs['units'] = "gm^-3"     
    p.close()
"""