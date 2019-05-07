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
from netCDF4 import Dataset
import matplotlib.pyplot as plt

start = time.time()

#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/GFDL/clt_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r') #Uni Laptop
dataset = Dataset('E:/University/University/MSc/Models/Data/GFDL/clt_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r') #Home PC

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
gfdl_tcc = np.take(gfdl_tcc, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
gfdl_tcc = np.mean(gfdl_tcc, axis=0) # average total cloud cover data over time
tcc = np.mean(gfdl_tcc, axis=1) # average total cloud cover data over longitude
gfdl_tcc_lat = np.vstack((lat, tcc)).T # Join the two lists as if they were two columns side by side, into a list of two elements each

#Select the months from 2010
gfdl_tciw = np.take(gfdl_tciw, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
gfdl_tciw = np.mean(gfdl_tciw, axis=0) # ice water content average over time
gfdl_tciw = np.mean(gfdl_tciw, axis=1) # ice water content average over longitude
#convert IWP to specific liquid water content (kg/kg) - divide by total air path
gfdl_tciw = gfdl_tciw / ap / 100

#Select the months from 2010
gfdl_tclw = np.take(gfdl_tclw, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
gfdl_tclw = np.mean(gfdl_tclw, axis=0) # liquid water content average over time
gfdl_tclw = np.mean(gfdl_tclw, axis=1) # liquid water content average over longitude
#convert IWP to specific liquid water content (kg/kg) - divide by total air path
gfdl_tclw = gfdl_tclw / ap / 100

# Join the two lists as if they were two columns side by side, into a list of two elements each
gfdl_tclw_lat = np.vstack((lat, gfdl_tclw)).T
gfdl_tciw_lat = np.vstack((lat, gfdl_tciw)).T
#----------------------------#

end = time.time()
print('Averaging data and creating combined arrays took:', end - start, 's')

#Select latitudes over the southern ocean
#gfdl_tcc_lat = gfdl_tcc_lat[gfdl_tcc_lat[:,0]>=-70]
#gfdl_tcc_lat = gfdl_tcc_lat[gfdl_tcc_lat[:,0]<=-50]

plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()

ax1.plot(gfdl_tcc_lat[:,0],gfdl_tcc_lat[:,1], '-r', label='Total Cloud Fraction')
ax2.plot(gfdl_tclw_lat[:,0],gfdl_tclw_lat[:,1], '-b', label='Specific Liquid Water Content (kg/kg)')
ax2.plot(gfdl_tciw_lat[:,0],gfdl_tciw_lat[:,1], '--b', label='Specific Liquid Water Content (kg/kg)')

#ax.axis('equal')
ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);
ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
          ncol=4, fancybox=True, shadow=True);
           
ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Specific Liquid and Ice Water Content (kg/kg)')

plt.title('Cloud Fraction and Specific Phase Content vs Latitude - GFDL.AM4 - 2010')

plt.grid(True)
plt.show()