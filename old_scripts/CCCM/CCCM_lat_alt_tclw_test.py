# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will create a datasset of specific cloud ice water content (kgm^-2) with latitude.
Data is stored in a 2D array cccm_tclw_lat 
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
from scipy import integrate
from sklearn.impute import SimpleImputer


input_path = 'D:/Downloads/CCCM/2006'
output_path = 'D:/Downloads/CCCM/tclw_output'

"""
p = h5py.File('../2006_CCCM_profile_variables.h5', 'r')
alt_labels = p['alt'][:]
p.close()
"""

os.chdir(input_path) #Uni laptop

#os.chdir('E:/University/University/MSc/Models/Data/CCCM/')  # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/Data/CCCM/Test') #Uni laptop
#os.chdir('//synthesis/E/University/University/MSc/Models/Data/CCCM/2006') #Uni laptop

def unique_lats( lat ):
    lat = np.array(lat)
    lat[:] = [(round(v*2,0)/2-90)*-1 for v in lat] #convert colatitude to latitude and raound to nearest 0.5 degrees
    return np.unique(lat)

def preprocess_tclw_data( tclw, lat ):
    
    tclw = np.array(tclw)
    #lat = np.array(lat)

    tclw  = np.transpose(tclw)
    tclw  = tclw[24:]
    tclw  = np.transpose(tclw)
          
    tclw[tclw > 20] = np.nan #Set the large 'fill values' in the data to nan before averaging  
    imp = SimpleImputer(missing_values=np.nan, strategy='mean')
    imp.fit(np.transpose(tclw))  
    a = imp.transform(np.transpose(tclw))
    tclw = np.transpose(a)
    
    lat = np.array(lat)
    lat[:] = [(round(v*2,0)/2-90)*-1 for v in lat] #convert colatitude to latitude and raound to nearest 0.5 degrees
    
    # a_tclw = tclw / 1000 #convert g to kg
    a_lat = np.vstack(lat)
    
    # Join the two lists
    combined = np.hstack((a_lat, tclw))
    
    # Add a column for every additional column, -1 will sort by the first column
    combined.sort(axis=0) #sort the combined array by latitude
    
    # total bin size of latitude entires empty list
    cccm_tclw_lat = np.empty((0,113),dtype=float)
    
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
            cccm_tclw_lat = np.concatenate((cccm_tclw_lat, average)) #stack new row
            # Reset the subtotal
            subtotal = np.zeros([0, 113])
            # Set the current latitude
            current_lat = item[0]
        
        # Add the next value to the subtotal
        subtotal=np.concatenate((subtotal, np.array([item[1:]]))) #stack the item on the bottom of the subtotal array
    
    # Catch the last entry in the for loop
    average = np.array([np.nanmean(subtotal, axis = 0)])
    cccm_tclw_lat = np.concatenate((cccm_tclw_lat, average))

    return cccm_tclw_lat


files = os.listdir()
file_groups = {}
# Group files into months
for file in files:
    # Example file name CER-NEWS_CCCM_Aqua-FM3-MODIS-CAL-CS_RelB1_905905.20060731.hdf
    # Extract the year month part of the file name
    group_name = "CCCM_tclw_lat_alt-" + file[-12:-6] + ".h5";
    # Add all the files of the same year month into the same list
    if group_name not in file_groups:
        file_groups[group_name] = []
    file_groups[group_name] += [file]

for output_file, input_files in file_groups.items():

    # Load every file in the directory
    tclw = np.empty((0,137), dtype=float) # create a blank array to add cloud ice water content data
    lat = [] # create a blank array to add lat values
    
    for input_file in input_files: 

        print(input_file)
        
        # Load the file
        f = SD.SD(input_file)
        # Get the cloud ice water content data as a list. (lat, alt) 'units': 'grams per cubic meter'
        tclw = np.concatenate((tclw,f.select('Liquid water content profile used').get()))
        lat = lat+f.select('Colatitude of subsatellite point at surface at observation').get().tolist()
    
    if len(lat) != tclw.shape[1]:
        exit('Invalid sizes of lat and tclw data')

    tclw_data = preprocess_tclw_data(tclw, lat)
    unique = unique_lats( lat )
    
    #print(test)
    #print(unique)


#    plt.subplots()
#    plt.contourf(lat1, cccm121_alt, np.transpose(cccm85_enhanced_lwc))
#    plt.show()


    
    with h5py.File(output_path + "/" + output_file, 'w') as p:
        p.create_dataset('lw_alt_lat', data=preprocess_tclw_data(tclw, lat))
        p.close()
    

#f.close()

