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

###############################################################################

cccm21_cloud_free_area = []
cccm34_phase = []
cccm52_cloud_fraction_profile = []
cccm67_cloudsat_lwc = []
cccm70_cloudsat_iwc = []
cccm85_enhanced_lwc = []
cccm86_enhanced_iwc = []

start = time.time()

# The directory where your HDF files are stored
os.chdir('d:/Downloads/CCCM')  # Home PC

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)

#    cccm21_cloud_free_area = cccm21_cloud_free_area + f.select('Cloud free area percent coverage (CALIPSO-CloudSat)').get().tolist()

#    cccm34_phase = cccm34_phase + f.select('Mean group cloud particle phase from MODIS radiance (3.7)').get().tolist()

#    cccm52_cloud_fraction_profile = cccm52_cloud_fraction_profile + f.select('Cloud fraction profile').get().tolist()

#    cccm67_cloudsat_lwc = cccm67_cloudsat_lwc + f.select('Mean CloudSat radar only liquid water content').get().tolist()

#    cccm85_enhanced_lwc = cccm85_enhanced_lwc + f.select('Liquid water content profile used').get().tolist()

    cccm86_enhanced_iwc = cccm86_enhanced_iwc + f.select('Ice water content profile used').get().tolist()

end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()

#cccm21_cloud_free_area = np.array(cccm21_cloud_free_area)
#cccm34_phase = np.array(cccm34_phase)
#cccm52_cloud_fraction_profile = np.array(cccm52_cloud_fraction_profile)
#cccm67_cloudsat_lwc = np.array(cccm67_cloudsat_lwc)
#cccm70_cloudsat_iwc = np.array(cccm70_cloudsat_iwc)
#cccm85_enhanced_lwc = np.array(cccm85_enhanced_lwc)
cccm86_enhanced_iwc = np.array(cccm86_enhanced_iwc)
end = time.time()
print('Create arrays took:', end - start, 's')



###############################################################################
import h5py
os.chdir('//synthesis/e/University/University/MSc/Models/climate-analysis/CCCM')
f = h5py.File('2010_raw_data.h5', 'r')
cccm21_cloud_free_area = f['cccm21_cloud_free_area'][:]
cccm34_phase = f['cccm34_phase'][:]
cccm52_cloud_fraction_profile = f['cccm52_cloud_fraction_profile'][:]
cccm67_cloudsat_lwc = f['cccm67_cloudsat_lwc'][:]
cccm70_cloudsat_iwc = f['cccm70_cloudsat_iwc'][:]

f.close()
# specify path and file name to create 
with h5py.File('2010_raw_data.h5', 'w') as p:
    p.create_dataset('cccm21_cloud_free_area', data=cccm21_cloud_free_area)
    p.create_dataset('cccm34_phase', data=cccm34_phase)
    p.create_dataset('cccm52_cloud_fraction_profile', data=cccm52_cloud_fraction_profile)
    p.create_dataset('cccm67_cloudsat_lwc', data=cccm67_cloudsat_lwc)
    p.create_dataset('cccm70_cloudsat_iwc', data=cccm70_cloudsat_iwc)

    p.close()

with h5py.File('2010_raw_cccm85_enhanced_lwc.h5', 'w') as p:

    p.create_dataset('cccm85_enhanced_lwc', data=cccm85_enhanced_lwc)

    p.close()
    

###############################################################################


os.chdir('e:/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets')  # Home PC
f = h5py.File('2010_CCCM_profile_variables.h5', 'r')

lat = f['lat'][:]
cccm123_alt = f['alt'][:]
cccm121_alt = f['alt_c'][:]
cccm124_alt = f['alt_t'][:]
f.close()

cccm86_enhanced_iwc[cccm86_enhanced_iwc > 20] = np.nan

###############################################################################

#fit all nan values to average

imp = SimpleImputer(missing_values=np.nan, strategy='mean')
imp.fit(np.transpose(cccm86_enhanced_iwc))  
a = imp.transform(np.transpose(cccm86_enhanced_iwc))
cccm86_enhanced_iwc = np.transpose(a)



###############################################################################

cccm86_enhanced_iwc_alt = np.nanmean(cccm86_enhanced_iwc, axis = 0)
cccm86_enhanced_iwc_lat = np.nanmean(cccm86_enhanced_iwc, axis = 1)

###############################################################################

# Join the two lists as if they were two columns side by side, into a list of two elements each
combined = np.vstack((lat, cccm86_enhanced_iwc_lat)).T
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

cccm86_enhanced_iwc_lat = cccm_tciw_lat
cccm86_enhanced_iwc_alt = np.vstack((cccm123_alt,cccm86_enhanced_iwc_alt)).T


os.chdir('//synthesis/e/University/University/MSc/Models/climate-analysis/CCCM')
with h5py.File('2010_cccm86_enhanced_iwc.h5', 'w') as p:

    p.create_dataset('cccm86_enhanced_iwc', data=cccm86_enhanced_iwc)
    p.create_dataset('cccm86_enhanced_iwc_alt', data=cccm86_enhanced_iwc_alt)
    p.create_dataset('cccm86_enhanced_iwc_lat', data=cccm86_enhanced_iwc_lat)

    p.close()    









###############################################################################
"""

os.chdir('e:/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets')  # Home PC
f = h5py.File('2010_CCCM_profile_variables.h5', 'r')


pressure_g = f['pressure_g_alt'][:]
pressure_so = f['pressure_so_alt'][:]
temp_g = f['temp_g_alt'][:]
temp_so = f['temp_so_alt'][:]
air_density_g = f['air_density_g'][:]
air_density_so = f['air_density_so'][:]

os.chdir('e:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')  # Home PC
f = h5py.File('2010_CCCM.h5', 'r')

cccm81b_cloud_area_enhanced =  f['tcc'][:]
cccm85_lwc_alt = f['lw'][:]
cccm86_iwc_alt = f['iw'][:]

cccm85_lwc_temp = f['lw_t'][:]
cccm86_iwc_temp = f['iw_t'][:]

cccm85_lwc_temp_so = f['lw_t_so'][:]
cccm86_iwc_temp_so = f['iw_t_so'][:]
"""    
    
    
    