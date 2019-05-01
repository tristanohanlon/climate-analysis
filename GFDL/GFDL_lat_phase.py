# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 18:40:37 2019

@author: Tristan O'Hanlon

This will create a dataset of time averaged global total cloud cover with latitude.
Time period is from 01.1980 - 12.2014
Data is stored in the 2D arrays: 

gfdl_tcc_lat
gfdl_tclw_lat
gfdl_tciw_lat

[:,0] = latitude
[:,1] = cloud fraction
"""
import time
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.pyplot as plt

start = time.time()

# Uni Laptop
#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/GFDL/clt_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
# Home PC
dataset = Dataset('E:/University/University/MSc/Models/Data/GFDL/clt_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')

lat = dataset.variables['lat'][:] #Extract latitude data
gfdl_tcc = dataset.variables['clt'][:] #Extract total cloud cover, keyed to time, lon and lat
        
dataset = Dataset('E:/University/University/MSc/Models/Data/GFDL/clivi_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
gfdl_tciw = dataset.variables['clivi'][:] #Extract ice water content (kg/m^2), keyed to time, lon and lat

dataset = Dataset('E:/University/University/MSc/Models/Data/GFDL/lwp_AERmon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
gfdl_tclw = dataset.variables['lwp'][:] #Extract liquid water content (kg/m^2), keyed to time, lon and lat

end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()

#Select the months from 2010
#gfdl_tcc = np.take(gfdl_tcc, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
gfdl_tcc = np.mean(gfdl_tcc, axis=0) # average total cloud cover data over time
tcc = np.mean(gfdl_tcc, axis=1) # average total cloud cover data over longitude
gfdl_tcc_lat = np.vstack((lat, tcc)).T # Join the two lists as if they were two columns side by side, into a list of two elements each

#Select the months from 2010
#gfdl_tciw = np.take(gfdl_tciw, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
gfdl_tciw = np.mean(gfdl_tciw, axis=0) # ice water content average over time
gfdl_tciw = np.mean(gfdl_tciw, axis=1) # ice water content average over longitude

#Select the months from 2010
#gfdl_tclw = np.take(gfdl_tciw, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
gfdl_tclw = np.mean(gfdl_tclw, axis=0) # liquid water content average over time
gfdl_tclw = np.mean(gfdl_tclw, axis=1) # liquid water content average over longitude

# Scale ice and liquid water content to total cloud cover
tclw_data = (gfdl_tclw/(gfdl_tciw+gfdl_tclw))*tcc
tciw_data = (gfdl_tciw/(gfdl_tciw+gfdl_tclw))*tcc

# Join the two lists as if they were two columns side by side, into a list of two elements each
gfdl_tclw_lat = np.vstack((lat, tclw_data)).T
gfdl_tciw_lat = np.vstack((lat, tciw_data)).T
#----------------------------#

end = time.time()
print('Averaging data and creating combined arrays took:', end - start, 's')

#Select latitudes over the southern ocean
#gfdl_tcc_lat = gfdl_tcc_lat[gfdl_tcc_lat[:,0]>=-70]
#gfdl_tcc_lat = gfdl_tcc_lat[gfdl_tcc_lat[:,0]<=-50]

plt.figure()
fig, ax = plt.subplots()
ax.plot(gfdl_tcc_lat[:,0],gfdl_tcc_lat[:,1], '-r', label='Total Cloud Fraction')
#ax.plot(gfdl_tclw_lat[:,0],gfdl_tclw_lat[:,1], '-b', label='Liquid Water Fraction')
#ax.plot(gfdl_tciw_lat[:,0],gfdl_tciw_lat[:,1], '--b', label='Ice Water Fraction')

#ax.axis('equal')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud and Phase Fraction')
plt.title('Total Cloud Fraction and Phase vs Latitude GFDL 1979 - Present')

plt.grid(True)
plt.show()