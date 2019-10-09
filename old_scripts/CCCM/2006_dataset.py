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
import matplotlib.pyplot as plt

###############################################################################

cccm21_cloud_free_area = []
cccm34_phase = []
cccm52_cloud_fraction_profile = []

start = time.time()

# The directory where your HDF files are stored
os.chdir('d:/Downloads/CCCM/2006')  # Home PC

# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)

    cccm21_cloud_free_area = cccm21_cloud_free_area + f.select('Cloud free area percent coverage (CALIPSO-CloudSat)').get().tolist()

    cccm34_phase = cccm34_phase + f.select('Mean group cloud particle phase from MODIS radiance (3.7)').get().tolist()

    cccm52_cloud_fraction_profile = cccm52_cloud_fraction_profile + f.select('Cloud fraction profile').get().tolist()

end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()

cccm21_cloud_free_area = np.array(cccm21_cloud_free_area)
cccm34_phase = np.array(cccm34_phase)
cccm52_cloud_fraction_profile = np.array(cccm52_cloud_fraction_profile)

end = time.time()
print('Create arrays took:', end - start, 's')



############################################################################### Save raw data

# specify path and file name to create 
os.chdir('e:/University/University/MSc/Models/Data/CCCM/raw_data')

with h5py.File('2006_raw_data.h5', 'w') as p:
    p.create_dataset('cccm21_cloud_free_area', data=cccm21_cloud_free_area)
    p.create_dataset('cccm34_phase', data=cccm34_phase)
    p.create_dataset('cccm52_cloud_fraction_profile', data=cccm52_cloud_fraction_profile)

    p.close()
    
############################################################################### Load raw data

"""
import h5py
os.chdir('//synthesis/e/University/University/MSc/Models/Data/CCCM/raw_data')
f = h5py.File('2006_raw_data.h5', 'r')
cccm21_cloud_free_area = f['cccm21_cloud_free_area'][:]
cccm34_phase = f['cccm34_phase'][:]
cccm52_cloud_fraction_profile = f['cccm52_cloud_fraction_profile'][:]

f.close()   

""" 

############################################################################### Get Southern Ocean Data - get lat first


# Join the two lists as if they were two columns side by side, into a list of two elements each
a = np.vstack(lat)
combined = np.hstack((a, cccm52_cloud_fraction_profile))
#print ("combined")
#print (combined)

# Add a column for every additional column, -1 will sort by the first column
combined = combined[np.lexsort(np.transpose(combined)[:-1])]
#print ("sorted")
#print (combined)

#Select latitudes over the southern ocean
combined = combined[combined[:,0]>=-70]
combined = combined[combined[:,0]<=-50]

#Split the combined array into just the lw data, eliminating the first coloumn of latitude
cccm52_cloud_fraction_profile_so = combined[:,1:114]
#alt = alt[0:137] #scale alt if necessary


############################################################################### Reduce and convert cloud free area to cloud area

cccm52_cloud_fraction_profile_alt = np.vstack((cccm121_alt, np.nanmean(cccm52_cloud_fraction_profile, axis = 0))).T
cccm52_cloud_fraction_profile_alt_so = np.vstack((cccm121_alt, np.nanmean(cccm52_cloud_fraction_profile_so, axis = 0))).T

cccm21_cloud_free_area = (100 - cccm21_cloud_free_area) / 100


############################################################################### Reduce cloud free area

# Join the two lists as if they were two columns side by side, into a list of two elements each
combined = np.vstack((lat, cccm21_cloud_free_area)).T
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
cccm_tcc_lat = np.empty((averages_total,2),dtype=float)

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
        cccm_tcc_lat[i] = [current_lat, average]
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
cccm_tcc_lat[i] = [current_lat, average]

cccm21_cloud_area_fraction = cccm_tcc_lat
os.chdir('e:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')  # Home PC
f = h5py.File('2006_CCCM.h5', 'r')

cccm81b_cloud_area_enhanced =  f['tcc'][:]
f.close()

############################################################################### Save reduced data


os.chdir('//synthesis/e/University/University/MSc/Models/Data/CCCM/raw_data')
with h5py.File('2006_cloud_fractions.h5', 'w') as p:

    p.create_dataset('cccm21_cloud_area_fraction', data=cccm21_cloud_area_fraction)
    p.create_dataset('cccm52_cloud_fraction_profile_alt', data=cccm52_cloud_fraction_profile_alt)
    p.create_dataset('cccm52_cloud_fraction_profile_alt_so', data=cccm52_cloud_fraction_profile_alt_so)
    p.create_dataset('cccm81b_cloud_area_enhanced', data=cccm81b_cloud_area_enhanced)

    p.close()
    
############################################################################### Load reduced data
  
"""    
os.chdir('//synthesis/e/University/University/MSc/Models/Data/CCCM/raw_data')
f = h5py.File('2006_cloud_fractions.h5', 'r')

cccm21_cloud_area_fraction = f['cccm21_cloud_area_fraction'][:]
cccm52_cloud_fraction_profile_alt = f['cccm52_cloud_fraction_profile_alt'][:]
cccm52_cloud_fraction_profile_alt_so = f['cccm52_cloud_fraction_profile_alt_so'][:]
#cccm81b_cloud_area_enhanced = f['cccm81b_cloud_area_enhanced'][:]
f.close()
"""
############################################################################### Load reduced phase data

os.chdir('e:/University/University/MSc/Models/Data/CCCM/raw_data')  # Home PC
f = h5py.File('2006_cccm85_enhanced_lwc.h5', 'r')
cccm85_enhanced_lwc_lat = f['cccm85_enhanced_lwc_lat'][:]
cccm85_enhanced_lwc_alt = f['cccm85_enhanced_lwc_alt'][:]
cccm85_enhanced_lwc_alt_so = f['cccm85_enhanced_lwc_alt_so'][:]

f.close()

os.chdir('e:/University/University/MSc/Models/Data/CCCM/raw_data')  # Home PC
f = h5py.File('2006_cccm86_enhanced_iwc.h5', 'r')
cccm86_enhanced_iwc_lat = f['cccm86_enhanced_iwc_lat'][:]
cccm86_enhanced_iwc_alt = f['cccm86_enhanced_iwc_alt'][:]
cccm86_enhanced_iwc_alt_so = f['cccm86_enhanced_iwc_alt_so'][:]

f.close()

############################################################################### Create phase fraction data


tclw_frac = np.vstack((cccm85_enhanced_lwc_lat[:,0], (cccm85_enhanced_lwc_lat[:,1] / (cccm85_enhanced_lwc_lat[:,1] + cccm86_enhanced_iwc_lat[:,1])) * cccm21_cloud_area_fraction[:,1])).T
tciw_frac = np.vstack((cccm86_enhanced_iwc_lat[:,0], (cccm86_enhanced_iwc_lat[:,1] / (cccm85_enhanced_lwc_lat[:,1] + cccm86_enhanced_iwc_lat[:,1])) * cccm21_cloud_area_fraction[:,1])).T

"""
fig, ax = plt.subplots()
ax.plot(cccm85_enhanced_lwc_lat[:,0], tclw_frac, '-k')
#ax.plot(cccm85_enhanced_lwc_lat[:,0], tclw_frac81, '-b')
ax.plot(cccm85_enhanced_lwc_lat[:,0], tciw_frac, '--k')
ax.plot(cccm21_cloud_free_area[:,0], cccm21_cloud_free_area[:,1], '--b')
"""

lw_frac = np.vstack((cccm85_enhanced_lwc_alt[:,0], (cccm85_enhanced_lwc_alt[:,1] / (cccm85_enhanced_lwc_alt[:,1] + cccm86_enhanced_iwc_alt[:,1])) * cccm52_cloud_fraction_profile_alt[:,1])).T
iw_frac = np.vstack((cccm86_enhanced_iwc_alt[:,0], (cccm86_enhanced_iwc_alt[:,1] / (cccm85_enhanced_lwc_alt[:,1] + cccm86_enhanced_iwc_alt[:,1])) * cccm52_cloud_fraction_profile_alt[:,1])).T

"""
fig, ax = plt.subplots()
ax.plot(lw_frac, cccm52_cloud_fraction_profile_alt[:,0], '-k')
#ax.plot(cccm85_enhanced_lwc_lat[:,0], tclw_frac81, '-b')
ax.plot(iw_frac, cccm52_cloud_fraction_profile_alt[:,0], '--k')
ax.plot(cccm52_cloud_fraction_profile_alt[:,1], cccm52_cloud_fraction_profile_alt[:,0], '--b')
"""

lw_frac_so = np.vstack((cccm85_enhanced_lwc_alt_so[:,0], (cccm85_enhanced_lwc_alt_so[:,1] / (cccm85_enhanced_lwc_alt_so[:,1] + cccm86_enhanced_iwc_alt_so[:,1])) * cccm52_cloud_fraction_profile_alt_so[:,1])).T
iw_frac_so = np.vstack((cccm86_enhanced_iwc_alt_so[:,0], (cccm86_enhanced_iwc_alt_so[:,1] / (cccm85_enhanced_lwc_alt_so[:,1] + cccm86_enhanced_iwc_alt_so[:,1])) * cccm52_cloud_fraction_profile_alt_so[:,1])).T



############################################################################### Load Previous data
"""
os.chdir('e:/University/University/MSc/Models/Data/CCCM/raw_data')  # Home PC
f = h5py.File('2006_lat.h5', 'r')

lat = f['lat'][:]
f.close()  
"""

os.chdir('e:/University/University/MSc/Models/climate-analysis/CCCM/raw_datasets')  # Home PC
f = h5py.File('2006_CCCM_profile_variables.h5', 'r')

cccm123_alt = f['alt'][:]
cccm121_alt = f['alt_c'][:]
cccm124_alt = f['alt_t'][:]
pressure_g = f['pressure_g_alt'][:]
pressure_so = f['pressure_so_alt'][:]
temp_g = f['temp_g_alt'][24:137]
temp_so = f['temp_so_alt'][24:137]
air_density_g = f['air_density_g'][:]
air_density_so = f['air_density_so'][:]
f.close()  
 
os.chdir('e:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')  # Home PC
f = h5py.File('2006_CCCM.h5', 'r')

tclw = f['tclw'][:]
tciw = f['tciw'][:]

cccm85_specific_lwc_alt = f['lw'][:]
cccm86_specific_iwc_alt = f['iw'][:]

cccm85_specific_lwc_alt_so = f['lw_so'][:]
cccm86_specific_iwc_alt_so = f['iw_so'][:]

cccm85_specific_lwc_temp = f['lw_t'][:]
cccm86_specific_iwc_temp = f['iw_t'][:]

cccm85_specific_lwc_temp_so = f['lw_t_so'][:]
cccm86_specific_iwc_temp_so = f['iw_t_so'][:]
f.close()   

############################################################################### Create temperature data


cf_frac_temp = np.vstack((temp_g[:, 1], cccm52_cloud_fraction_profile_alt[:,1])).T
lw_frac_temp = np.vstack((temp_g[:, 1], lw_frac[:,1])).T
iw_frac_temp = np.vstack((temp_g[:, 1], iw_frac[:,1])).T

cf_frac_temp_so = np.vstack((temp_so[:, 1], cccm52_cloud_fraction_profile_alt_so[:,1])).T
lw_frac_temp_so = np.vstack((temp_so[:, 1], lw_frac_so[:,1])).T
iw_frac_temp_so = np.vstack((temp_so[:, 1], iw_frac_so[:,1])).T

############################################################################### Create new datasets
  
os.chdir('e:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')  # Home PC
with h5py.File('2006_data_new.h5', 'w') as p:

    p.create_dataset('lat', data=lat)
    p.create_dataset('alt', data=cccm121_alt)
    p.create_dataset('air_density_g', data=air_density_g)
    p.create_dataset('air_density_so', data=air_density_so)
  
    p.create_dataset('tcc', data=cccm21_cloud_area_fraction)
    
    p.create_dataset('tclw', data=tclw)
    p.create_dataset('tciw', data=tciw)
    
    p.create_dataset('tclw_gcm3', data=cccm85_enhanced_lwc_lat)
    p.create_dataset('tciw_gcm3', data=cccm86_enhanced_iwc_lat)
   
    p.create_dataset('tclw_frac', data=tclw_frac)
    p.create_dataset('tciw_frac', data=tciw_frac)
    
    p.create_dataset('cf', data=cccm52_cloud_fraction_profile_alt)
    p.create_dataset('cf_so', data=cccm52_cloud_fraction_profile_alt_so)
    
    p.create_dataset('lw_frac', data=lw_frac)
    p.create_dataset('lw_frac_so', data=lw_frac_so)
    p.create_dataset('iw_frac', data=iw_frac)
    p.create_dataset('iw_frac_so', data=iw_frac_so)
    
    p.create_dataset('lw', data=cccm85_specific_lwc_alt)
    p.create_dataset('lw_so', data=cccm85_specific_lwc_alt_so)
    p.create_dataset('iw', data=cccm86_specific_iwc_alt)
    p.create_dataset('iw_so', data=cccm86_specific_iwc_alt_so)

    p.create_dataset('temp', data=temp_g)
    p.create_dataset('temp_so', data=temp_so)
    p.create_dataset('pressure', data=pressure_g)
    p.create_dataset('pressure_so', data=pressure_so)
   
    p.create_dataset('lw_t', data=cccm85_specific_lwc_temp)
    p.create_dataset('lw_t_so', data=cccm85_specific_lwc_temp_so)
    p.create_dataset('iw_t', data=cccm86_specific_iwc_temp)
    p.create_dataset('iw_t_so', data=cccm86_specific_iwc_temp_so)
    
    p.create_dataset('cf_t', data=cf_frac_temp)
    p.create_dataset('cf_t_so', data=cf_frac_temp_so)
    p.create_dataset('lw_frac_temp', data=lw_frac_temp)
    p.create_dataset('lw_frac_temp_so', data=lw_frac_temp_so)
    p.create_dataset('iw_frac_temp', data=iw_frac_temp)
    p.create_dataset('iw_frac_temp_so', data=iw_frac_temp_so)

    p.close()


############################################################################### Test plots
  
"""
fig, ax = plt.subplots()
ax.plot(tclw_frac[:,0], tclw_frac[:, 1], '-k')
ax.plot(tciw_frac[:,0], tciw_frac[:, 1], '--k')
ax.plot(cccm21_cloud_free_area[:,0], cccm21_cloud_free_area[:,1], '--b')


fig, ax = plt.subplots()
ax.plot(lw_frac_so[:,1], lw_frac_so[:, 0], '-k')
ax.plot(iw_frac_so[:,1], iw_frac_so[:, 0], '--k')
ax.plot(cccm52_cloud_fraction_profile_alt_so[:,1], cccm52_cloud_fraction_profile_alt_so[:,0], '-b')


fig, ax = plt.subplots()
ax.plot(lw_frac_temp[:,0], lw_frac_temp[:, 1], '-k')
ax.plot(iw_frac_temp[:,0], iw_frac_temp[:, 1], '--k')
ax.plot(cf_frac_temp[:,0], cf_frac_temp[:,1], '-b')

fig, ax = plt.subplots()
ax.plot(lw_frac_temp[:,0], lw_frac_temp[:, 1], '-k')
ax.plot(lw_frac_temp_so[:,0], lw_frac_temp_so[:, 1], '--k')

"""    