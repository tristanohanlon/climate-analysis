# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon


"""
import time
import numpy as np
import os
from pyhdf import SD
import h5py

#Empty lists
lat =[]
alt = []
alt_c = []
alt_t = []
temp = []
pressure = []
cff = []

###############################################################################

start = time.time()

# The directory where your HDF files are stored
#os.chdir('C:/Users/toha006/University/University/MSc/Models/Data/CCCM/2010')  # Uni Laptop
os.chdir('E:/University/University/MSc/Models/Data/CCCM/2006')  # Home PC

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    # Get the altitude levels (m) which corespond to the ice water content data. Only store the last file's values
    alt = f.select('Irradiance layer center height profile').get().tolist() #137 levels for ice and liquid water profiles
    # Get the altitude levels (km) which corespond to the cloud fraction data. Only store the last file's values, already in km.
    alt_c = f.select('Layer center height profile (clouds and aerosol)').get().tolist() #113 levels for clouds and aerosols
    # Get the altitude levels (m) which corresponds to the pressure and temperature data. Only store the last file's values
    alt_t = f.select('Irradiance level height profile').get().tolist() #138 altitude levels pressure. 
    # Get the latitude data as a list
    lat = lat+f.select('Colatitude of subsatellite point at surface at observation').get().tolist()
    #Get the pressure (hPa) data which is a (lat, alt_t) array per file. Axis = 0 is added to for each file.
    pressure = pressure+f.select('Pressure profile').get().tolist()
    #Get the temperature (K) data which is a (lat, alt_t) array per file. Axis = 0 is added to for each file.
    temp = temp+f.select('Temperature profile').get().tolist()
    cff = cff+f.select('Cloud free area percent coverage (CALIPSO-CloudSat)').get().tolist() #25536 values each referenced to latitude and longitude

    end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()

lat[:] = [(round(v*2,0)/2-90)*-1 for v in lat]
#print("round lats")
lat = np.array(lat)
alt = np.array(alt) / 1000 # Convert the altitude list to a numpy array and convery m to km
alt_t = np.array(alt_t) / 1000 # Convert the altitude list to a numpy array and convery m to km
alt_c = np.array(alt_c)
pressure = np.array(pressure)
temp = np.array(temp)
cff = np.array(cff) 

#Set the large 'fill values' in the data to nan before averaging
cff[cff > 100] = None        
pressure[pressure > 1100] = None   
temp[temp > 400] = None   

# Average the ice water content over latitude and longitude for each altitude level 
pressure_g = np.nanmean(pressure, axis=0)
temp_g = np.nanmean(temp, axis=0)

temp_g_alt = np.vstack((alt_t, temp_g)).T
pressure_g_alt = np.vstack((alt_t, pressure_g)).T

end = time.time()
print('Create arrays and averaging data sets creation took:', end - start, 's')


###############################################################################

start = time.time()

#Calculate global density profile and divide IWC to get specific IWC at each profile level

air_density_g = [] #create empty list

#calculate air density (gm^-3) at each altitude layer
air_density_g = (pressure_g * 100) / (286.9 * temp_g) * 1000
air_density_g = air_density_g[1:138]

end = time.time()
print('Create global air density data set took:', end - start, 's')

###############################################################################

start = time.time()

#---Southern Ocean Temperature Profile---#

# Join the two lists as if they were two columns side by side, into a list of two elements each
a = np.vstack(lat)
combined_t = np.hstack((a, temp))
#print ("combined")
#print (combined)

# Add a column for every additional column, -1 will sort by the first column
combined_t = combined_t[np.lexsort(np.transpose(combined_t)[:-1])]
#print ("sorted")
#print (combined)

#Select latitudes over the southern ocean
combined_t = combined_t[combined_t[:,0]>=-70]
combined_t = combined_t[combined_t[:,0]<=-50]

#Split the combined array into just the iw data, eliminating the first coloumn of latitude
temp_so = combined_t[:,1:139]
#alt_t = alt_t[0:137] #scale alt if necessary

temp_so = np.nanmean(temp_so, axis=0)
temp_so_alt = np.vstack((alt_t, temp_so)).T


#---Southern Ocean Pressure Profile---#

# Join the two lists as if they were two columns side by side, into a list of two elements each
combined_p = np.hstack((a, pressure))
#print ("combined")
#print (combined)

# Add a column for every additional column, -1 will sort by the first column
combined_p = combined_p[np.lexsort(np.transpose(combined_p)[:-1])]
#print ("sorted")
#print (combined)

#Select latitudes over the southern ocean
combined_p = combined_p[combined_p[:,0]>=-70]
combined_p = combined_p[combined_p[:,0]<=-50]

#Split the combined array into just the iw data, eliminating the first coloumn of latitude
pressure_so = combined_p[:,1:139]
#alt_t = alt_t[0:137] #scale alt if necessary

pressure_so = np.nanmean(pressure_so, axis=0)
pressure_so_alt = np.vstack((alt_t, pressure_so)).T

end = time.time()
print('Create southern ocean data set took:', end - start, 's')

###############################################################################

start = time.time()

#Calculate southern ocean density profile and divide IWC to get specific IWC at each profile level

air_density_so = [] #create empty list

#calculate air density (gm^-3) at each altitude layer
air_density_so = (pressure_so * 100) / (286.9 * temp_so) * 1000
air_density_so = air_density_so[1:138]

end = time.time()
print('Create southern ocean air density data set took:', end - start, 's')


###############################################################################

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
cffb = 1 - cffa
cff = tciw[:,1]*cffb

#Select latitudes over the southern ocean
co = cff_lat[cff_lat[:,0]>=-70]
co = co[co[:,0]<=-50]

#Split the combined array into just the iw data, eliminating the first coloumn of latitude
cff_so = co[:,1]
cff_so = 1 - cff_so

###############################################################################
import h5py

os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets')
# specify path and file name to create 
with h5py.File('2006_CCCM_profile_variables.h5', 'w') as p:
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    