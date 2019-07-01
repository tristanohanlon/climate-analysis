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
tcc = dataset.variables['clt'][:] #Extract total cloud cover, keyed to time, lon and lat
        
dataset = Dataset('E:/University/University/MSc/Models/Data/GFDL/clivi_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
tciw = dataset.variables['clivi'][:] #Extract ice water content (kg/m^2), keyed to time, lon and lat

dataset = Dataset('E:/University/University/MSc/Models/Data/GFDL/lwp_AERmon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
tclw = dataset.variables['lwp'][:] #Extract liquid water content (kg/m^2), keyed to time, lon and lat

end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()


#gfdl_tcc = np.take(gfdl_tcc, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0) #Select the months from 2010

#Data from July 2006 to April 2011 corresponding to CCCM Data
tcc = np.take(tcc, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 
tcc = np.mean(tcc, axis=0) # average total cloud cover data over time
tcc = np.mean(tcc, axis=1) # average total cloud cover data over longitude
tcc = np.vstack((lat, tcc)).T # Join the two lists as if they were two columns side by side, into a list of two elements each


#gfdl_tciw = np.take(gfdl_tciw, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0) #Select the months from 2010

#Data from July 2006 to April 2011 corresponding to CCCM Data
tciw = np.take(tciw, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 
tciw = np.mean(tciw, axis=0) # ice water content average over time
tciw = np.mean(tciw, axis=1) # ice water content average over longitude


#gfdl_tclw = np.take(gfdl_tclw, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)#Select the months from 2010

#Data from July 2006 to April 2011 corresponding to CCCM Data
tclw = np.take(tclw, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 
tclw = np.mean(tclw, axis=0) # liquid water content average over time
tclw = np.mean(tclw, axis=1) # liquid water content average over longitude


# Join the two lists as if they were two columns side by side, into a list of two elements each
tclw = np.vstack((lat, tclw)).T
tciw = np.vstack((lat, tciw)).T
#----------------------------#

end = time.time()
print('Averaging data and creating combined arrays took:', end - start, 's')

#Select latitudes over the southern ocean
#gfdl_tcc_lat = gfdl_tcc_lat[gfdl_tcc_lat[:,0]>=-70]
#gfdl_tcc_lat = gfdl_tcc_lat[gfdl_tcc_lat[:,0]<=-50]

plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()

ax1.plot(tcc[:,0],tcc[:,1], '-r', label='Total Cloud Fraction')
ax2.plot(tclw[:,0],tclw[:,1], '-b', label='Liquid Water Content')
ax2.plot(tciw[:,0],tciw[:,1], '--b', label='Ice Water Content')

#ax.axis('equal')
ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);
ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
          ncol=4, fancybox=True, shadow=True);
           
ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Liquid and Ice Water Content ($kgm^{-2}$)')

plt.title('Cloud Fraction and Phase Content vs Latitude - GFDL.AM4 - July 2006 to April 2011')

plt.grid(True)
plt.show()