# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

This will create a dataset of time averaged global cloud cover, liquid water cloud fraction and ice water cloud fraction with altitude.
The code can select either global or southern ocean data.
Data is stored in the 2D arrays: 

gfdl_tcc_alt
gfdl_tclw_alt
gfdl_tciw_alt

[:,0] = alt
[:,1] = cloud fraction

The data has already been scaled and offsetted.
"""
import time
import numpy as np
import math as math
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import h5py
import os
from scipy import integrate

start = time.time()
os.chdir('E:/University/University/MSc/Models/Data/GFDL/') #Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/Data/GFDL/') #Uni Laptop

dataset = Dataset('clw_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
lw = dataset.variables['clw'][:] #Extract Mass Fraction of Cloud Liquid Water (kg/kg), keyed to time, level, lon and lat (420, 33, 180, 288)
lat = dataset.variables['lat'][:] #Extract latitude, (180,)

dataset = Dataset('cli_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
iw = dataset.variables['cli'][:] #Extract Mass Fraction of Cloud Ice (kg/kg), keyed to time, level, lon and lat (420, 33, 180, 288)

end = time.time()
print('Importing data took:', end - start, 's')

#----------------------------#

#Select the months from 07.2006 to 04.2011
lw = np.take(lw, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 
iw = np.take(iw, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 
#Average over time
lw = np.mean(lw, axis=0)
iw = np.mean(iw, axis=0)
#Average over longitude
lw = np.mean(lw, axis=-1)
iw = np.mean(iw, axis=-1)

#Select Southern ocean Latitudes

lw_so = np.vstack((lat, lw)).T #creates a (180,34) array
lw_so = lw_so[lw_so[:,0]>=-70]
lw_so = lw_so[lw_so[:,0]<=-50]
lw_so = lw_so[:,1:34] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude

iw_so = np.vstack((lat, iw)).T #creates a (180,34) array
iw_so = iw_so[iw_so[:,0]>=-70]
iw_so = iw_so[iw_so[:,0]<=-50]
iw_so = iw_so[:,1:34] #Split the combined array into just the tciw data, eliminating the first coloumn of latitude

lw_alt_lat = lw
iw_alt_lat = iw

lw_alt_lat_so = np.transpose(lw_so)
iw_alt_lat_so = np.transpose(iw_so)

###############################################################################

#Import old data
os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL/reduced_datasets')
b = h5py.File('07.2006_04.2011_gfdl.h5', 'r')

tcc = b['tcc'][:] # 0-1
tclw = b['tclw'][:] #kgm^-2
tciw = b['tciw'][:] #kgm^-2

#---gfdl Global Profile---#

gfdl_tcc_alt = b['cf'][:] # 0-1
gfdl_tclw_alt = b['lw'][:] #kg/kg
gfdl_tciw_alt = b['iw'][:] #kg/kg
gfdl_temp_alt = b['temp'][:] #K
gfdl_pressure_alt = b['pressure'][:] #hPa

gfdl_tcc_temp = b['cf_t'][:] # 0-1
gfdl_tclw_temp = b['lw_t'][:] #kg/kg
gfdl_tciw_temp = b['iw_t'][:] #kg/kg


#---gfdl Southern Ocean Profile---#

gfdl_tcc_alt_so = b['cf_so'][:] # 0-1
gfdl_tclw_alt_so = b['lw_so'][:] #kg/kg
gfdl_tciw_alt_so = b['iw_so'][:] #kg/kg
gfdl_temp_alt_so = b['temp_so'][:] #K
gfdl_pressure_alt_so = b['pressure_so'][:] #hPa

gfdl_tcc_temp_so = b['cf_t_so'][:] # 0-1
gfdl_tclw_temp_so = b['lw_t_so'][:] #kg/kg
gfdl_tciw_temp_so = b['iw_t_so'][:] #kg/kg

alt = gfdl_tcc_alt[:,0]
lat_so = lat[20:40]



###############################################################################

os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL/reduced_datasets') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GFDL/reduced_datasets') #Uni Laptop

with h5py.File('07.2006_04.2011_gfdla.h5', 'w') as p:

    p.create_dataset('lat', data=lat)
    p.create_dataset('lat_so', data=lat_so)  
    p.create_dataset('alt', data=alt_fixed)
    
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('tclw', data=tclw)
    p.create_dataset('tciw', data=tciw)  
    
    p.create_dataset('cf', data=gfdl_tcc_alt)
    p.create_dataset('cf_so', data=gfdl_tcc_alt_so)
    p.create_dataset('lw', data=gfdl_tclw_alt)
    p.create_dataset('lw_so', data=gfdl_tclw_alt_so)
    p.create_dataset('iw', data=gfdl_tciw_alt)
    p.create_dataset('iw_so', data=gfdl_tciw_alt_so)

    p.create_dataset('temp', data=gfdl_temp_alt)
    p.create_dataset('temp_so', data=gfdl_temp_alt_so)
    p.create_dataset('pressure', data=gfdl_pressure_alt)
    p.create_dataset('pressure_so', data=gfdl_pressure_alt_so)
    
    p.create_dataset('lw_alt_lat', data=a)
    p.create_dataset('iw_alt_lat', data=b)
    p.create_dataset('lw_alt_lat_so', data=c)
    p.create_dataset('iw_alt_lat_so', data=d)
    
    p.create_dataset('cf_t', data=gfdl_tcc_temp)
    p.create_dataset('cf_t_so', data=gfdl_tcc_temp_so)
    p.create_dataset('lw_t', data=gfdl_tclw_temp)
    p.create_dataset('lw_t_so', data=gfdl_tclw_temp_so)
    p.create_dataset('iw_t', data=gfdl_tciw_temp)
    p.create_dataset('iw_t_so', data=gfdl_tciw_temp_so)

    p.close()


fig, ax = plt.subplots()
ax.contourf(lat_so, alt_fixed[0:13], c[0:13])
ax.colorbar()















