# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

This will create a dataset of time averaged global total cloud cover, liquid water cloud fraction and ice water cloud fraction with latitude.
Data is stored in the 2D arrays: 

ecmwf_tcc_lat
ecmwf_tclw_lat
ecmwf_tciw_lat

[:,0] = latitude
[:,1] = cloud fraction
"""
import time
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.pyplot as plt

start = time.time()
dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/1979-03.2019_ECMWF_amon_tcc_tciw_tclw.nc', 'r')
#lon = dataset.variables['longitude'][:] #Extract longitude data
lat = dataset.variables['latitude'][:] #Extract latitude data
tcc = dataset.variables['tcc'][:] #Extract total cloud cover, keyed to time, lon and lat
tciw = dataset.variables['tciw'][:] #Extract ice water content (kg/m^2), keyed to time, lon and lat
tclw = dataset.variables['tclw'][:] #Extract liquid water content (kg/m^2), keyed to time, lon and lat
end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()
#--------total cloud cover------------#

#Select the months from 2010
#tcc = np.take(tcc, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)

tc = np.mean(tcc, axis=0) # total cloud cover averaged over time
atcc = np.mean(tc, axis=1) # total cloud cover averaged over longitude
ecmwf_tcc_lat = np.vstack((lat, atcc)).T # Join the two lists as if they were two columns side by side, into a list of two elements each

#--------ice water content--------#

tciw = np.take(tciw, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
tci = np.mean(tciw, axis=0) # ice water content average over time
atciw = np.mean(tci, axis=1) # ice water content average over longitude

#--------liquid water content--------#
tclw = np.take(tclw, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
tcl = np.mean(tclw, axis=0) # liquid water content average over time
atclw = np.mean(tcl, axis=1) # liquid water content average over longitude

# Scale ice and liquid water content to total cloud cover
tclw_data = (atclw / (atciw+atclw) ) * atcc
tciw_data = (atciw / (atciw+atclw) ) * atcc

# Join the two lists as if they were two columns side by side, into a list of two elements each
ecmwf_tclw_lat = np.vstack((lat, tclw_data)).T
ecmwf_tciw_lat = np.vstack((lat, tciw_data)).T
#----------------------------#
end = time.time()
print('Averaging data and creating combined arrays took:', end - start, 's')

#Select latitudes over the southern ocean
#ecmwf_tcc_lat = ecmwf_tcc_lat[ecmwf_tcc_lat[:,0]>=-70]
#ecmwf_tcc_lat = ecmwf_tcc_lat[ecmwf_tcc_lat[:,0]<=-50]

plt.figure()
fig, ax = plt.subplots()
ax.plot(ecmwf_tcc_lat[:,0],ecmwf_tcc_lat[:,1], '-r', label='Cloud Fraction')
ax.plot(ecmwf_tclw_lat[:,0],ecmwf_tclw_lat[:,1], '-b', label='Liquid Water Fraction')
ax.plot(ecmwf_tciw_lat[:,0],ecmwf_tciw_lat[:,1], '--b', label='Ice Water Fraction')

#ax.axis('equal')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud and Phase Fraction')
plt.title('Cloud Fraction and Phase vs Latitude ECMWF 2010')

plt.grid(True)
plt.show()