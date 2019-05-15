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
from netCDF4 import Dataset
import h5py
import os
import matplotlib.pyplot as plt

start = time.time()

# Uni Laptop
#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/1979-03.2019_ECMWF_amon_tcc_tciw_tclw.nc', 'r')
# Home PC
dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF/1979-03.2019_ECMWF_amon_tcc_tciw_tclw.nc', 'r')

#lon = dataset.variables['longitude'][:] #Extract longitude data
lat = dataset.variables['latitude'][:] #Extract latitude data
tcc = dataset.variables['tcc'][:] #Extract total cloud cover, keyed to time, lon and lat
tciw = dataset.variables['tciw'][:] #Extract ice water path (kg/m^2), keyed to time, lon and lat
tclw = dataset.variables['tclw'][:] #Extract liquid water path (kg/m^2), keyed to time, lon and lat

end = time.time()
print('Importing data from files to lists took:', end - start, 's')

start = time.time()
#--------total cloud cover------------#
#tcc = np.take(tcc, [319, 320, 321, 322, 323, 324], axis=0) #Select the months July to December from 2006
#tcc = np.take(tcc, [325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336], axis=0) #Select the months from 2007
#tcc = np.take(tcc, [337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348], axis=0) #Select the months from 2008
#tcc = np.take(tcc, [349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360], axis=0) #Select the months from 2009
#tcc = np.take(tcc, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0) #Select the months from 2010
#tcc = np.take(tcc, [374,375,376], axis=0) #Select the months Feb to April from 2011

#Data from July 2006 to April 2011 corresponding to CCCM Data
tcc = np.take(tcc, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 

ecmwf_tcc_lat = np.mean(tcc, axis=0) # total cloud cover averaged over time
ecmwf_tcc_lat = np.mean(ecmwf_tcc_lat, axis=1) # total cloud cover averaged over longitude
ecmwf_tcc_lat = np.vstack((lat, ecmwf_tcc_lat)).T # Join the two lists as if they were two columns side by side, into a list of two elements each

#--------ice water content--------#

#tciw = np.take(tciw, [319, 320, 321, 322, 323, 324], axis=0) #Select the months July to December from 2006
#tciw = np.take(tciw, [325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336], axis=0) #Select the months from 2007
#tciw = np.take(tciw, [337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348], axis=0) #Select the months from 2008
#tciw = np.take(tciw, [349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360], axis=0) #Select the months from 2009
#tciw = np.take(tciw, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0) #Select the months from 2010
#tciw = np.take(tciw, [374,375,376], axis=0) #Select the months Feb to April from 2011

#Data from July 2006 to April 2011 corresponding to CCCM Data
tciw = np.take(tciw, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 

ecmwf_tciw_lat = np.mean(tciw, axis=0) # ice water content average over time
ecmwf_tciw_lat= np.mean(ecmwf_tciw_lat , axis=1) # ice water content average over longitude


#--------liquid water content--------#

#tclw = np.take(tclw, [319, 320, 321, 322, 323, 324], axis=0) #Select the months July to December from 2006
#tclw = np.take(tclw, [325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336], axis=0) #Select the months from 2007
#tclw = np.take(tclw, [337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348], axis=0) #Select the months from 2008
#tclw = np.take(tclw, [349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360], axis=0) #Select the months from 2009
#tclw = np.take(tclw, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0) #Select the months from 2010
#tclw = np.take(tclw, [374,375,376], axis=0) #Select the months Feb to April from 2011

#Data from July 2006 to April 2011 corresponding to CCCM Data
tclw = np.take(tclw, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 

ecmwf_tclw_lat = np.mean(tclw, axis=0) # liquid water content average over time
ecmwf_tclw_lat = np.mean(ecmwf_tclw_lat, axis=1) # liquid water content average over longitude


# Join the two lists as if they were two columns side by side, into a list of two elements each
ecmwf_tclw_lat = np.vstack((lat, ecmwf_tclw_lat)).T
ecmwf_tciw_lat = np.vstack((lat, ecmwf_tciw_lat)).T
#----------------------------#
end = time.time()
print('Averaging data and creating combined arrays took:', end - start, 's')

plt.figure()
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.plot(ecmwf_tcc_lat[:,0],ecmwf_tcc_lat[:,1], '-r', label='Cloud Fraction')
ax2.plot(ecmwf_tclw_lat[:,0],ecmwf_tclw_lat[:,1], '-b', label='Liquid Water Content ($kgm^{-2}$)')
ax2.plot(ecmwf_tciw_lat[:,0],ecmwf_tciw_lat[:,1], '--b', label='Ice Water Content ($kgm^{-2}$)')

#ax.axis('equal')
ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
          ncol=4, fancybox=True, shadow=True);
          
ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Phase Content ($kgm^{-2}$)')

plt.title('Cloud Fraction and Phase Content vs Latitude ECMWF 2006 - 2011')

plt.grid(True)
plt.show()

###############################################################################

os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets') #Uni Laptop

with h5py.File('2006_2011_ECMWF_global_latitude.h5', 'w') as p:
    
    p.create_dataset('tcc', data=ecmwf_tcc_lat)
    p.create_dataset('tclw', data=ecmwf_tclw_lat)
    p.create_dataset('tciw', data=ecmwf_tciw_lat)
    
    p.close()
