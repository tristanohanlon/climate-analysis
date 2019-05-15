# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

This will create a dataset of time averaged global cloud cover, liquid water cloud fraction and ice water cloud fraction with altitude.
The code can select either global or southern ocean data.

Data is stored in the 2D arrays: 

ecmwf_tcc_plevel
ecmwf_tclw_plevel
ecmwf_tciw_plevel
ecmwf_temp_plevel

[:,0] = alt
[:,1] = cloud fraction

The data has already been scaled and offsetted.
"""
import time
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

#https://www.mide.com/pages/air-pressure-at-altitude-calculator
#https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
alt = [39.92, 38.01, 35.4, 31.99, 30.09, 28.50, 26.06, 23.90, 20.64, 18.50, 16.23, 14.81, 13.64, 12.66, 11.531, 11.026, 8.024, 7.291, 6.652, 6.049, 5.467, 4.901, 4.351, 3.817, 3.299, 2.798, 2.314, 2.078, 1.845, 1.615, 1.389, 1.167, 0.949, 0.734, 0.522, 0.314, 0.108]
alt = np.array(alt)

# Uni Laptop
#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/pressure_levels/2010_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
# Home PC
dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/pressure_levels/2010_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')

start = time.time()

lat = dataset.variables['latitude'][:] #Extract latitude data (721)
plevel = dataset.variables['level'][:] #Extract pressure (millibars) level data (37 levels)
tcc = dataset.variables['cc'][:] #Extract fraction of cloud cover, keyed to time, lat and lon (12, 37, 721, 1440)
tciw = dataset.variables['ciwc'][:] #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
tclw = dataset.variables['clwc'][:] #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
temp = dataset.variables['t'][:] #Extract air temperature profile content (K), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

end = time.time()
print('Importing data took:', end - start, 's')

#---------------tcc-------------#
start = time.time()

tcc = np.array(tcc)
tcc = np.mean(tcc, axis=0) #Average fraction of cloud cover over time
tcc = np.mean(tcc, axis=-1) #Average fraction of cloud cover over longitude

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
tcc_so = np.vstack((lat, tcc)).T #creates a (721,38) array

#Select latitudes over the southern ocean
tcc_so = tcc_so[tcc_so[:,0]>=-70]
tcc_so = tcc_so[tcc_so[:,0]<=-50]

#Split the combined array into just the tcc data, eliminating the first coloumn of latitude
tcc_so = tcc_so[:,1:38]

tcc = np.mean(tcc, axis=-1) #Average fraction of cloud cover over latitude
tcc_so = np.mean(tcc_so, axis=0) #Average southern ocean fraction of cloud cover over latitude

end = time.time()
print('Averaging tcc data took:', end - start, 's')

#---------------tciw-------------#

start = time.time()

tciw = np.array(tciw)
tciw = np.mean(tciw, axis=0) #Average specific cloud ice water content (kg/kg) over time
tciw = np.mean(tciw, axis=-1) #Average specific cloud ice water content (kg/kg) over longitude

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
tciw_so = np.vstack((lat, tciw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
tciw_so = tciw_so[tciw_so[:,0]>=-70]
tciw_so = tciw_so[tciw_so[:,0]<=-50]

#Split the combined array into just the tciw data, eliminating the first coloumn of latitude
tciw_so = tciw_so[:,1:38]

tciw = np.mean(tciw, axis=-1) #Average specific cloud ice water content (kg/kg) over latitude
tciw_so = np.mean(tciw_so, axis=0) #Average southern ocean specific cloud ice water content (kg/kg) over latitude

end = time.time()
print('Averaging tciw data took:', end - start, 's')

#---------------tclw-------------#

start = time.time()

tclw = np.array(tclw)
tclw = np.mean(tclw, axis=0) #Average specific cloud liquid water content (kg/kg) over time
tclw = np.mean(tclw, axis=-1) #Average specific cloud liquid water content (kg/kg) over longitude

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
tclw_so = np.vstack((lat, tclw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
tclw_so = tclw_so[tclw_so[:,0]>=-70]
tclw_so = tclw_so[tclw_so[:,0]<=-50]

#Split the combined array into just the tclw data, eliminating the first coloumn of latitude
tclw_so = tclw_so[:,1:38]

tclw = np.mean(tclw, axis=-1) #Average specific cloud liquid water content (kg/kg) over latitude
tclw_so = np.mean(tclw_so, axis=0) #Average southern ocean specific cloud liquid water content (kg/kg) over latitude

end = time.time()
print('Averaging tclw data took:', end - start, 's')

#---------------temperature-------------#

start = time.time()

temp = np.array(temp)
temp = np.mean(temp, axis=0) #Average air temperature over time
temp = np.mean(temp, axis=-1) #Average air temperature over longitude

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
temp_so = np.vstack((lat, temp)).T #creates a (721,38) array

#Select latitudes over the southern ocean
temp_so = temp_so[temp_so[:,0]>=-70]
temp_so = temp_so[temp_so[:,0]<=-50]

#Split the combined array into just the temp data, eliminating the first coloumn of latitude
temp_so = temp_so[:,1:38]

temp = np.mean(temp, axis=-1) #Average air temperature over latitude
temp_so = np.mean(temp_so, axis=0)  #Average southern ocean air temperature over latitude

end = time.time()
print('Averaging temp data took:', end - start, 's')

#---------------combine-------------#

start = time.time()

# Join the two lists as if they were two columns side by side, into a list of two elements each
ecmwf_tcc_alt = np.vstack((alt, tcc)).T 
ecmwf_tclw_alt = np.vstack((alt, tclw)).T
ecmwf_tciw_alt = np.vstack((alt, tciw)).T
ecmwf_temp_alt = np.vstack((alt, temp)).T
ecmwf_plevel_alt = np.vstack((alt, plevel)).T

ecmwf_tcc_temp = np.vstack((temp, tcc)).T 
ecmwf_tclw_temp = np.vstack((temp, tclw)).T
ecmwf_tciw_temp = np.vstack((temp, tciw)).T

ecmwf_tcc_alt_so = np.vstack((alt, tcc_so)).T 
ecmwf_tclw_alt_so = np.vstack((alt, tclw_so)).T
ecmwf_tciw_alt_so = np.vstack((alt, tciw_so)).T
ecmwf_temp_alt_so = np.vstack((alt, temp_so)).T
ecmwf_plevel_alt_so = np.vstack((alt, plevel)).T

ecmwf_tcc_temp_so = np.vstack((temp_so, tcc_so)).T 
ecmwf_tclw_temp_so = np.vstack((temp_so, tclw_so)).T
ecmwf_tciw_temp_so = np.vstack((temp_so, tciw_so)).T

#ecmwf_tcc_plevel = np.vstack((plevel, tcc)).T 
#ecmwf_tclw_plevel = np.vstack((plevel, tclw)).T
#ecmwf_tciw_plevel = np.vstack((plevel, tciw)).T

end = time.time()
print('Creating the combined arrays took:', end - start, 's')

#----------------------------#

plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()
ax1.plot(ecmwf_tcc_alt[:,1],ecmwf_tcc_alt[:,0], '-r', label='Fraction Cloud Cover')
ax2.plot(ecmwf_tclw_alt[:,1]*10000,ecmwf_tclw_alt[:,0], '-b', label='Specific Cloud Liquid Water Content')
ax2.plot(ecmwf_tciw_alt[:,1]*10000,ecmwf_tciw_alt[:,0], '--b', label='Specific Cloud Ice Water Content')

#ax.axis('equal')
ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);
#ax.axis('equal')
ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
          ncol=4, fancybox=True, shadow=True);
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Specific Cloud Liquid and Ice Water Content (kg/kg) x $10^{-4}$')
ax1.set_ylabel('Altitude (km)')
plt.title('Southern Ocean Cloud Fraction and Phase vs Altitude ECMWF 2010')
#plt.gca().invert_yaxis()

plt.grid(True)
plt.show()

###############################################################################
import h5py
import os
os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets') #Uni Laptop

with h5py.File('2010_ECMWF_global_profile.h5', 'w') as p:
    
    p.create_dataset('cf', data=ecmwf_tcc_alt)
    p.create_dataset('cf_so', data=ecmwf_tcc_alt_so)
    p.create_dataset('lw', data=ecmwf_tclw_alt)
    p.create_dataset('lw_so', data=ecmwf_tclw_alt_so)
    p.create_dataset('iw', data=ecmwf_tciw_alt)
    p.create_dataset('iw_so', data=ecmwf_tciw_alt_so)

    p.create_dataset('temp', data=ecmwf_temp_alt)
    p.create_dataset('temp_so', data=ecmwf_temp_alt_so)
    p.create_dataset('pressure', data=ecmwf_plevel_alt)
    p.create_dataset('pressure_so', data=ecmwf_plevel_alt)

    p.create_dataset('cf_t', data=ecmwf_tcc_temp)
    p.create_dataset('cf_t_so', data=ecmwf_tcc_temp_so)
    p.create_dataset('lw_t', data=ecmwf_tclw_temp)
    p.create_dataset('lw_t_so', data=ecmwf_tclw_temp_so)
    p.create_dataset('iw_t', data=ecmwf_tciw_temp)
    p.create_dataset('iw_t_so', data=ecmwf_tciw_temp_so)

    p.close()

