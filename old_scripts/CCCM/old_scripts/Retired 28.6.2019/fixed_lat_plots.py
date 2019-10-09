# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will extract the global ice water content profile against altitude. 
The ice water content is averaged over all longitude and latitude at each altitude layer.
"""
import time
import numpy as np
import os
from pyhdf import SD
import h5py


#Import altitude and latitude data from reduced datasets

#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Uni Laptop
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets') #Home PC
p = h5py.File('2011_CCCM_profile_variables.h5', 'r')

alt = p['alt'][:]
lat = p['lat'][:]
air_density_g = p['air_density_g'][:]
air_density_so = p['air_density_so'][:]
alt_c = p['alt_c'][:]
alt_t = p['alt_t'][:]
temp_g_alt = p['temp_g_alt'][:]
pressure_g_alt = p['pressure_g_alt'][:]
temp_so_alt = p['temp_so_alt'][:]
pressure_so_alt = p['pressure_so_alt'][:]


os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Home PC
f = h5py.File('2011_CCCM.h5', 'r')

tcc = f['tcc'][:]
tclw = f['tclw'][:]
tciw = f['tciw'][:]
temp = f['temp'][:]
pressure = f['pressure'][:]

temp_so = f['temp_so'][:]
pressure_so = f['pressure_so'][:]

cf = f['cf'][:]
lw = f['lw'][:]
iw = f['iw'][:]

cf_so = f['cf_so'][:]
lw_so = f['lw_so'][:]
iw_so = f['iw_so'][:]

cf_t = f['cf_t'][:]
lw_t = f['lw_t'][:]
iw_t = f['iw_t'][:]

cf_t_so = f['cf_t_so'][:]
lw_t_so = f['lw_t_so'][:]
iw_t_so = f['iw_t_so'][:]


###############################################################################

cff = []
# The directory where your HDF files are stored
os.chdir('E:/University/University/MSc/Models/Data/CCCM/2011')  # Home PC

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the ice water content (gm^-3) data which is a (lat, alt) array per file. Axis = 0 is added to for each file.  
    cff = cff+f.select('Cloud free area percent coverage (CALIPSO-CloudSat)').get().tolist() # 137 levels, 25536 values each referenced to latitude and longitude

cff = np.array(cff) # Convert ice water content list to a numpy array

#Set the large 'fill values' in the data to nan before averaging
cff[cff > 100] = None        

# Join the two lists as if they were two columns side by side, into a list of two elements each
#a = np.vstack(lat)
combined = np.vstack((lat, cff)).T
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
cff_lat = np.empty((averages_total,2),dtype=float)

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
        #average = subtotal / number
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
        cff_lat[i] = [current_lat, average]
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
#average = subtotal / number
average = subtotal / number / 100
cff_lat[i] = [current_lat, average]

cffa = cff_lat[:,1]

cff = 1 - cffa

#Select latitudes over the southern ocean
co = cff_lat[cff_lat[:,0]>=-70]
co = co[co[:,0]<=-50]

#Split the combined array into just the iw data, eliminating the first coloumn of latitude
cff_so = co[:,1]
cff_so = 1 - cff_so

tclw = tclw[:,1]*cff
tclw = np.vstack((unique,tclw)).T
tciw = tciw[:,1]*cff
tciw = np.vstack((unique,tciw)).T


###############################################################################
import h5py

os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets')
# specify path and file name to create 
with h5py.File('2011_CCCM_profile_variablesa.h5', 'w') as p:
    p.create_dataset('lat', data=lat)
    p.create_dataset('alt', data=alt)
    p.create_dataset('alt_c', data=alt_c)
    p.create_dataset('alt_t', data=alt_t)
    p.create_dataset('temp_g_alt', data=temp_g_alt)
    p.create_dataset('pressure_g_alt', data=pressure_g_alt)
    p.create_dataset('temp_so_alt', data=temp_so_alt)
    p.create_dataset('pressure_so_alt', data=pressure_so_alt)
    p.create_dataset('air_density_g', data=air_density_g)
    p.create_dataset('air_density_so', data=air_density_so)  
    p.create_dataset('cff', data=cff) 
    p.create_dataset('cff_so', data=cff)  
    p.close()
    

os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') # Home PC

with h5py.File('2011_CCCMa.h5', 'w') as p:
    p.create_dataset('lat', data=lat)
    
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('tclw', data=tclw)
    p.create_dataset('tciw', data=tciw)
    
    p.create_dataset('alt', data=alt)
    p.create_dataset('alt_c', data=alt_c)
    p.create_dataset('alt_t', data=alt_t)

    p.create_dataset('air_density_g', data=air_density_g)
    p.create_dataset('air_density_so', data=air_density_so)

    p.create_dataset('cf', data=cf)
    p.create_dataset('cf_so', data=cf_so)
    p.create_dataset('lw', data=lw)
    p.create_dataset('lw_so', data=lw_so)
    p.create_dataset('iw', data=iw)
    p.create_dataset('iw_so', data=iw_so)

    p.create_dataset('temp', data=temp)
    p.create_dataset('temp_so', data=temp_so)
    p.create_dataset('pressure', data=pressure)
    p.create_dataset('pressure_so', data=pressure_so)

    p.create_dataset('cf_t', data=cf_t)
    p.create_dataset('cf_t_so', data=cf_t_so)
    p.create_dataset('lw_t', data=lw_t)
    p.create_dataset('lw_t_so', data=lw_t_so)
    p.create_dataset('iw_t', data=iw_t)
    p.create_dataset('iw_t_so', data=iw_t_so)

    p.close()    
