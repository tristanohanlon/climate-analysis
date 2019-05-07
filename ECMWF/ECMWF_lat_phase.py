# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

This will create a dataset of time averaged global total cloud cover, specific liquid water and specific ice water with latitude.
Data is stored in the 2D arrays: 

ecmwf_tcc_lat
ecmwf_tclw_lat
ecmwf_tciw_lat

[:,0] = latitude
[:,1] = cloud fraction or specific phase amount
"""
import time
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.pyplot as plt

start = time.time()

# Uni Laptop
dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/1979-03.2019_ECMWF_amon_tcc_tciw_tclw.nc', 'r')
# Home PC
#dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/1979-03.2019_ECMWF_amon_tcc_tciw_tclw.nc', 'r')

#lon = dataset.variables['longitude'][:] #Extract longitude data
lat = dataset.variables['latitude'][:] #Extract latitude data
#tcc = dataset.variables['tcc'][:] #Extract total cloud cover, keyed to time, lon and lat
tciw = dataset.variables['tciw'][:] #Extract ice water path (kg/m^2), keyed to time, lon and lat
tclw = dataset.variables['tclw'][:] #Extract liquid water path (kg/m^2), keyed to time, lon and lat

end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()
#--------total cloud cover------------#
"""
#Select the months from 2010
tcc = np.take(tcc, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
ecmwf_tcc_lat = np.mean(tcc, axis=0) # total cloud cover averaged over time
ecmwf_tcc_lat = np.mean(ecmwf_tcc_lat, axis=1) # total cloud cover averaged over longitude
ecmwf_tcc_lat = np.vstack((lat, ecmwf_tcc_lat)).T # Join the two lists as if they were two columns side by side, into a list of two elements each
"""
#--------ice water content--------#
#Select the months from 2010
tciw = np.take(tciw, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
ecmwf_tciw_lat = np.mean(tciw, axis=0) # ice water content average over time
ecmwf_tciw_lat= np.mean(ecmwf_tciw_lat , axis=1) # ice water content average over longitude
ecmwf_tciw_lat = ecmwf_tciw_lat / -ap #convert IWP to specific ice water content (kg/kg) by dividing by the total air path


#--------liquid water content--------#
#Select the months from 2010
tclw = np.take(tclw, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
ecmwf_tclw_lat = np.mean(tclw, axis=0) # liquid water content average over time
ecmwf_tclw_lat = np.mean(ecmwf_tclw_lat, axis=1) # liquid water content average over longitude
ecmwf_tclw_lat = ecmwf_tclw_lat / -ap #convert LWP to specific liquid water content (kg/kg)


# Join the two lists as if they were two columns side by side, into a list of two elements each
ecmwf_tclw_lat = np.vstack((lat, ecmwf_tclw_lat)).T
ecmwf_tciw_lat = np.vstack((lat, ecmwf_tciw_lat)).T
#----------------------------#
end = time.time()
print('Averaging data and creating combined arrays took:', end - start, 's')

#Select latitudes over the southern ocean
#ecmwf_tcc_lat = ecmwf_tcc_lat[ecmwf_tcc_lat[:,0]>=-70]
#ecmwf_tcc_lat = ecmwf_tcc_lat[ecmwf_tcc_lat[:,0]<=-50]

plt.figure()
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.plot(ecmwf_tcc_lat[:,0],ecmwf_tcc_lat[:,1], '-r', label='Cloud Fraction')
ax2.plot(ecmwf_tclw_lat[:,0],ecmwf_tclw_lat[:,1]*10000, '-b', label='Specific liquid Water Content (kg/kg)x$10^{-4}$')
ax2.plot(ecmwf_tciw_lat[:,0],ecmwf_tciw_lat[:,1]*10000, '--b', label='Specific Ice Water Content (kg/kg)x$10^{-4}$')

#ax.axis('equal')
ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
          ncol=4, fancybox=True, shadow=True);
          
ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Specific Phase Content (kg/kg)x$10^{-4}$')

plt.title('Cloud Fraction and Specific Phase Content vs Latitude ECMWF 2010')

plt.grid(True)
plt.show()