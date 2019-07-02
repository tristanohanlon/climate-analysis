# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon
"""

#---Import Latitude Variables---#

import h5py
import os
import numpy as np
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Uni Laptop
#os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Laptop

a = h5py.File('2007_2008_cccm.h5', 'r')
tcc = a['tcc'][:]
tclw = a['tclw'][:]
tciw = a['tciw'][:]
temp = a['temp'][:]
pressure = a['pressure'][:]

temp_so = a['temp_so'][:]
pressure_so = a['pressure_so'][:]

cf = a['cf'][:]
lw = a['lw'][:]
iw = a['iw'][:]

cf_so = a['cf_so'][:]
lw_so = a['lw_so'][:]
iw_so = a['iw_so'][:]

cf_t = a['cf_t'][:]
lw_t = a['lw_t'][:]
iw_t = a['iw_t'][:]

cf_t_so = a['cf_t_so'][:]
lw_t_so = a['lw_t_so'][:]
iw_t_so = a['iw_t_so'][:]

tclw_gcm3 = a['tclw_gcm3'][:]
tciw_gcm3 = a['tciw_gcm3'][:]
   
tclw_frac = a['tclw_frac'][:]
tciw_frac = a['tciw_frac'][:]
    
lw_frac = a['lw_frac'][:]
lw_frac_so = a['lw_frac_so'][:]
iw_frac = a['iw_frac'][:]
iw_frac_so = a['iw_frac_so'][:]
    
lw_frac_temp = a['lw_frac_temp'][:]
lw_frac_temp_so = a['lw_frac_temp_so'][:]
iw_frac_temp = a['iw_frac_temp'][:]
iw_frac_temp_so = a['iw_frac_temp_so'][:]

a.close()

os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets/old_data') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets/old_data') #Uni Laptop
a = h5py.File('2007_CCCM.h5', 'r')
tcc_enhanced1 = a['tcc'][:]
a.close()

a = h5py.File('2008_CCCM.h5', 'r')
tcc_enhanced2 = a['tcc'][:]
a.close()

############################################################################### tcc

combined = np.vstack((tcc_enhanced1, tcc_enhanced2))

#print("get unique lats")
unique = np.unique(tcc_enhanced1[:,0])
#print(unique)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]

# Averages of (lat, cloud ice water content) empty array
averages_total = unique.size
tcc_enhanced = np.empty((averages_total,2),dtype=float)

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
        current_lat = item[0];
    
    # If the lat is not the same as last time, then perform the average calc and reset everything
    if item[0] != current_lat:
        
        # Find the average value.
        average = subtotal / number
        # Append the average
        tcc_enhanced[i] = [current_lat, average]
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
tcc_enhanced[i] = [current_lat, average]

############################################################################### tclw


tclw_frac_enhanced = tclw_frac / tcc
tciw_frac_enhanced = tciw_frac / tcc

tclw_frac_enhanced = tclw_frac_enhanced * tcc_enhanced
tciw_frac_enhanced = tciw_frac_enhanced * tcc_enhanced

###############################################################################
    

import h5py
import os
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Uni Laptop
#os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Laptop

with h5py.File('2007_2008_cccm.h5', 'w') as p:
    
    p.create_dataset('tcc_enhanced', data=tcc_enhanced)
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('tclw', data=tclw)
    p.create_dataset('tciw', data=tciw)
    
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
    
    p.create_dataset('tclw_gcm3', data=tclw_gcm3)
    p.create_dataset('tciw_gcm3', data=tciw_gcm3)
   
    p.create_dataset('tclw_frac_enhanced', data=tclw_frac_enhanced)
    p.create_dataset('tciw_frac_enhanced', data=tciw_frac_enhanced)
    
    p.create_dataset('tclw_frac', data=tclw_frac)
    p.create_dataset('tciw_frac', data=tciw_frac)
    
    p.create_dataset('lw_frac', data=lw_frac)
    p.create_dataset('lw_frac_so', data=lw_frac_so)
    p.create_dataset('iw_frac', data=iw_frac)
    p.create_dataset('iw_frac_so', data=iw_frac_so)
 
    p.create_dataset('lw_frac_temp', data=lw_frac_temp)
    p.create_dataset('lw_frac_temp_so', data=lw_frac_temp_so)
    p.create_dataset('iw_frac_temp', data=iw_frac_temp)
    p.create_dataset('iw_frac_temp_so', data=iw_frac_temp_so)
    
    p.close()






















