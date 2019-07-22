# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will create a datasset of specific cloud ice water content (kgm^-2) with latitude.
Data is stored in a 2D array cccm_tcc_lat 
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
import gc



input_path = 'e:/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets/tcc_alt_lat'
output_path = 'e:/University/University/MSc/Models/climate-analysis/CCCM/test_data'


os.chdir(input_path) #Uni laptop

#os.chdir('E:/University/University/MSc/Models/Data/CCCM/')  # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/Data/CCCM/Test') #Uni laptop
#os.chdir('//synthesis/E/University/University/MSc/Models/Data/CCCM/2006') #Uni laptop


def preprocess_tcc_data( tcc, lat ):
    
    tcc = np.array(tcc)
    #lat = np.array(lat)
    
    a_lat = np.vstack(lat)
    
    # Join the two lists
    combined = np.hstack((a_lat, tcc))
    
    # Add a column for every additional column, -1 will sort by the first column
    combined.sort(axis=0) #sort the combined array by latitude
    
    # total bin size of latitude entires empty list
    cccm_tcc_lat = np.empty((0,113),dtype=float)
    
    # Current subtotal of current lat
    subtotal = np.zeros([0, 113])
    # Set the current lat to false
    current_lat = None
    
    # Iterate through all of the (lat, cloud ice water content) elements and subtotal the same lat values
    for item in combined:
    
        if current_lat is None:
            current_lat = item[0];
        
        # If the lat is not the same as last time, then perform the average calc and reset everything
        if item[0] != current_lat:
            
            # Find the average value.
            average = np.array([np.nanmean(subtotal, axis = 0)]) # exclude first row with zeros
            # Append the average
            cccm_tcc_lat = np.concatenate((cccm_tcc_lat, average)) #stack new row
            # Reset the subtotal
            subtotal = np.zeros([0, 113])
            # Set the current latitude
            current_lat = item[0]
        
        # Add the next value to the subtotal
        subtotal=np.concatenate((subtotal, np.array([item[1:]]))) #stack the item on the bottom of the subtotal array
    
    # Catch the last entry in the for loop
    average = np.array([np.nanmean(subtotal, axis = 0)])
    cccm_tcc_lat = np.concatenate((cccm_tcc_lat, average))

    return cccm_tcc_lat


files = os.listdir()
file_groups = {}
# Group files into months
for file in files:
    # Example file name CER-NEWS_CCCM_Aqua-FM3-MODIS-CAL-CS_RelB1_905905.20060731.hdf
    # Extract the year month part of the file name
    group_name = "2007-2008-tcc_lat_alt-" + ".h5";
    # Add all the files of the same year month into the same list
    if group_name not in file_groups:
        file_groups[group_name] = []
    file_groups[group_name] += [file]

for output_file, input_files in file_groups.items():

    # Load every file in the directory
    tcc = np.empty((0,113), dtype=float) # create a blank array to add cloud ice water content data
    lat = np.empty(0, dtype=float)
    
    for input_file in input_files: 
        os.chdir(input_path) #Uni laptop

        print(input_file)
        
        # Load the file
        f = h5py.File(input_file, 'r')
        # Get the cloud ice water content data as a list. (lat, alt) 'units': 'grams per cubic meter'
        tcc = np.concatenate((tcc,f['cf_alt_lat'][:]))
        f.close()
        
        os.chdir('e:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')  # Home PC
        f = h5py.File('2007_2008_cccm.h5', 'r')

        lat = np.concatenate((lat, f['tcc'][:,0]))

        f.close()

    
    tcc_data = preprocess_tcc_data(tcc, lat)
    
    
    with h5py.File(output_path + "/" + output_file, 'w') as p:
        p.create_dataset('cf_alt_lat', data=preprocess_tcc_data(tcc, lat))
        p.close()
    


"""
#---Import Data---#

#os.chdir('c:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')  # Uni     
os.chdir('e:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')  # Home PC
f = h5py.File('2007_2008_cccm.h5', 'r')

lat = f['tcc'][:,0]
alt = f['cf'][:,0]

f.close()

#os.chdir('c:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')  # Uni     
os.chdir('e:/University/University/MSc/Models/climate-analysis/CCCM/test_data')  # Home PC
f = h5py.File('2007-2008-tcc_lat_alt-.h5', 'r')

cf_alt_lat = f['cf_alt_lat'][:]

f.close()


#os.chdir('c:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')  # Uni     
os.chdir('e:/University/University/MSc/Models/climate-analysis/CCCM/test_data')  # Home PC
f = h5py.File('2007-2008-tciw_lat_alt-.h5', 'r')

iw_alt_lat = f['iw_alt_lat'][:]

f.close()

#os.chdir('c:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')  # Uni     
os.chdir('e:/University/University/MSc/Models/climate-analysis/CCCM/test_data')  # Home PC
f = h5py.File('2007-2008-tclw_lat_alt-.h5', 'r')

lw_alt_lat = f['lw_alt_lat'][:]

f.close()


#---Fit lw and iw---#

iw_alt_lat = np.transpose(iw_alt_lat)
iw_alt_lat = iw_alt_lat[:-24]
iw_alt_lat = np.transpose(iw_alt_lat)

lw_alt_lat = np.transpose(lw_alt_lat)
lw_alt_lat = lw_alt_lat[:-24]
lw_alt_lat = np.transpose(lw_alt_lat)


#---Fit nan values---#

imp = SimpleImputer(missing_values=np.nan, strategy='mean')
imp.fit(np.transpose(cccm85_enhanced_lwc))  
a = imp.transform(np.transpose(cccm85_enhanced_lwc))
cccm85_enhanced_lwc = np.transpose(a)


#---Create Fractions---#

lw_frac_alt_lat = (lw_alt_lat / (lw_alt_lat + iw_alt_lat)) * (cf_alt_lat)
iw_frac_alt_lat = (iw_alt_lat / (lw_alt_lat + iw_alt_lat)) * (cf_alt_lat)


#---Test plots---#


plt.subplots()
plt.contourf(lat, alt, np.transpose(cf_alt_lat))
plt.clim(0, 0.5)
plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.show()

plt.subplots()
plt.contourf(lat, alt, lw_frac_alt_lat)
plt.clim(0, 0.5)
plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.show()

plt.subplots()
plt.contourf(lat, alt, iw_frac_alt_lat)
plt.clim(0, 0.5)
plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.show()
"""