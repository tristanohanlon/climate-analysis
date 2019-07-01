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

dataset = Dataset('cl_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
cf = dataset.variables['cl'][:] #Extract cloud_area_fraction_in_atmosphere_layer, keyed to time, lev, lon and lat (420, 33, 180, 288)
lat = dataset.variables['lat'][:] #Extract latitude, (180,)

dataset = Dataset('clw_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
lw = dataset.variables['clw'][:] #Extract Mass Fraction of Cloud Liquid Water (kg/kg), keyed to time, level, lon and lat (420, 33, 180, 288)

dataset = Dataset('cli_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
iw = dataset.variables['cli'][:] #Extract Mass Fraction of Cloud Ice (kg/kg), keyed to time, level, lon and lat (420, 33, 180, 288)

dataset = Dataset('ta_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
temp = dataset.variables['ta'][:] #Extract air temperature, keyed to time, plev, lon and lat (420, 19, 180, 288)
plev = dataset.variables['plev'][:] #Extract air pressure corresponding to temperature (19,)

dataset = Dataset('pfull_AERmon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
pressure = dataset.variables['pfull'][:] #Extract air pressure, keyed to time, level, lon and lat (420, 33, 180, 288)

end = time.time()
print('Importing data took:', end - start, 's')

#----------------------------#

#Select the months from 2010
cf = np.take(cf, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 
lw = np.take(lw, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 
iw = np.take(iw, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 
temp = np.take(temp, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 
pressure = np.take(pressure, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 
#Average over time
cf = np.mean(cf, axis=0)
lw = np.mean(lw, axis=0)
iw = np.mean(iw, axis=0)
temp = np.mean(temp, axis=0)
pressure = np.mean(pressure, axis=0)
#Average over longitude
cf = np.mean(cf, axis=-1)
lw = np.mean(lw, axis=-1)
iw = np.mean(iw, axis=-1)
temp = np.mean(temp, axis=-1)
pressure = np.mean(pressure, axis=-1)

#Select Southern ocean Latitudes

cf_so = np.vstack((lat, cf)).T #creates a (180,34) array
cf_so = cf_so[cf_so[:,0]>=-70]
cf_so = cf_so[cf_so[:,0]<=-50]
cf_so = cf_so[:,1:34] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude

lw_so = np.vstack((lat, lw)).T #creates a (180,34) array
lw_so = lw_so[lw_so[:,0]>=-70]
lw_so = lw_so[lw_so[:,0]<=-50]
lw_so = lw_so[:,1:34] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude

iw_so = np.vstack((lat, iw)).T #creates a (180,34) array
iw_so = iw_so[iw_so[:,0]>=-70]
iw_so = iw_so[iw_so[:,0]<=-50]
iw_so = iw_so[:,1:34] #Split the combined array into just the tciw data, eliminating the first coloumn of latitude

temp_so = np.vstack((lat, temp)).T #creates a (180,34) array
temp_so = temp_so[temp_so[:,0]>=-70]
temp_so = temp_so[temp_so[:,0]<=-50]
temp_so = temp_so[:,1:34] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude

pressure_so = np.vstack((lat, pressure)).T #creates a (180,34) array
pressure_so = pressure_so[pressure_so[:,0]>=-70]
pressure_so = pressure_so[pressure_so[:,0]<=-50]
pressure_so = pressure_so[:,1:34] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude

#Average over latitude - change axis to 0 if getting southern ocean data
cf = np.mean(cf, axis=-1) / 100
lw = np.mean(lw, axis=-1)
iw = np.mean(iw, axis=-1)
temp = np.mean(temp, axis=-1)
pressure = np.mean(pressure, axis=-1) / 100

cf_so = np.mean(cf_so, axis=0) / 100
lw_so = np.mean(lw_so, axis=0)
iw_so = np.mean(iw_so, axis=0)
temp_so = np.mean(temp_so, axis=0)
pressure_so = np.mean(pressure_so, axis=0) / 100

###############################################################################

#---convert pressure levels to altitude---#

#https://www.mide.com/pages/air-pressure-at-altitude-calculator
#https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
alt_t = np.empty((pressure.size,1),dtype=float)
alt_t_so = np.empty((pressure.size,1),dtype=float)
alt_p = np.empty((pressure.size,1),dtype=float)
alt_ts = np.empty((pressure.size,1),dtype=float)
alt_ts_so = np.empty((pressure.size,1),dtype=float)


# Iterate through all of the temp elements (troposphere h < 11km)
i = 0
for item in pressure:
    newalt = (288.19 - 288.08*((item/1012.9)**(1/5.256)))/6.49
    alt_t[i] = [newalt]
    i+=1


# Iterate through all of the temp elements (troposphere h < 11km)
i = 0
for item in pressure_so:
    newalt = (288.19 - 288.08*((item/1012.9)**(1/5.256)))/6.49
    alt_t_so[i] = [newalt]
    i+=1


# Iterate through all of the pressure elements (lower stratosphere 11km < h <25km)
i = 0
for item in pressure:
    newalt = (1.73 - math.log(item/226.50))/0.157
    alt_p[i] = [newalt]
    i+=1


# Iterate through all of the temp elements (upper stratosphere  h > 25km)
i = 0
for item in pressure:
    newalt = (216.6*((item/24.88)**(1/-11.388)) - 141.94)/2.99
    alt_ts[i] = [newalt]
    i+=1


# Iterate through all of the temp elements (upper stratosphere  h > 25km)
i = 0
for item in pressure_so:
    newalt = (216.6*((item/24.88)**(1/-11.388)) - 141.94)/2.99
    alt_ts_so[i] = [newalt]
    i+=1


#manually adjust alt and alt_so arrays usinf alt_p and alt_ts
alt = alt_t
alt_so = alt_t_so 
    

#---create temperature profile---#
#troposphere from pressure data
#stratosphere assumed from Earth Atmosphere Model 
#https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html

temp_t = 288.19 - 6.49*alt
temp_s = 141.94 + 2.99*alt
#manually add on 216.69K values to altitudes above 11km to 25km 

temp = temp_t
temp_so = temp_t

###############################################################################
"""
#---------------combine-------------#

#alt_so = np.transpose(alt_so) #may need these to adjust combinations
#alt = np.transpose(alt)
# Join the two lists as if they were two columns side by side, into a list of two elements each
gfdl_tcc_alt = np.vstack((alt, cf)).T
gfdl_tclw_alt = np.vstack((alt, lw)).T
gfdl_tciw_alt = np.vstack((alt, iw)).T
gfdl_temp_alt = np.hstack((alt, temp))
gfdl_pressure_alt = np.vstack((alt, pressure)).T

gfdl_tcc_temp = np.vstack((temp, cf)).T
gfdl_tclw_temp = np.vstack((temp, lw)).T
gfdl_tciw_temp = np.vstack((temp, iw)).T

gfdl_tcc_alt_so = np.vstack((alt_so, cf_so)).T 
gfdl_tclw_alt_so = np.vstack((alt_so, lw_so)).T
gfdl_tciw_alt_so = np.vstack((alt_so, iw_so)).T
gfdl_temp_alt_so = np.hstack((alt_so, temp_so))
gfdl_pressure_alt_so = np.vstack((alt_so, pressure_so)).T

gfdl_tcc_temp_so = np.vstack((temp_so, cf_so)).T 
gfdl_tclw_temp_so = np.vstack((temp_so, lw_so)).T
gfdl_tciw_temp_so = np.vstack((temp_so, iw_so)).T

#----------------------------#


###############################################################################
import h5py
import os
os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL/reduced_datasets') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GFDL/reduced_datasets') #Uni Laptop

with h5py.File('07.2006_04.2011_gfdl.h5', 'w') as p:
    
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

    p.create_dataset('cf_t', data=gfdl_tcc_temp)
    p.create_dataset('cf_t_so', data=gfdl_tcc_temp_so)
    p.create_dataset('lw_t', data=gfdl_tclw_temp)
    p.create_dataset('lw_t_so', data=gfdl_tclw_temp_so)
    p.create_dataset('iw_t', data=gfdl_tciw_temp)
    p.create_dataset('iw_t_so', data=gfdl_tciw_temp_so)

    p.close()
"""


















