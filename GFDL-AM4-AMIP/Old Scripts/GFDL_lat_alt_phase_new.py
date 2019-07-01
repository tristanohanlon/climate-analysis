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
os.chdir('E:/University/University/MSc/Models/Data/CMIP6/gfdl_am4') #Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/Data/GFDL/') #Uni Laptop
dataset = Dataset('cl_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
cf = dataset.variables['cl'][:] #Extract Cloud Fraction, keyed to time, level, lon and lat (420, 33, 180, 288)

dataset = Dataset('clw_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
lw = dataset.variables['clw'][:] #Extract Mass Fraction of Cloud Liquid Water (kg/kg), keyed to time, level, lon and lat (420, 33, 180, 288)

dataset = Dataset('cli_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
iw = dataset.variables['cli'][:] #Extract Mass Fraction of Cloud Ice (kg/kg), keyed to time, level, lon and lat (420, 33, 180, 288)

dataset = Dataset('ta_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
temp = dataset.variables['ta'][:] #Extract Temperature, keyed to time, level, lon and lat (420, 19, 180, 288)

end = time.time()
print('Importing data took:', end - start, 's')

#----------------------------#

#Select the months from 07.2006 to 04.2011
cf = np.take(cf, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 
lw = np.take(lw, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 
iw = np.take(iw, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 
temp = np.take(temp, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 

#Average over time
cf = np.mean(cf, axis=0)
lw = np.mean(lw, axis=0)
iw = np.mean(iw, axis=0)
temp = np.mean(temp, axis=0)
#Average over longitude
cf_alt_lat = np.mean(cf, axis=-1)
lw_alt_lat = np.mean(lw, axis=-1)
iw_alt_lat = np.mean(iw, axis=-1)
temp_alt_lat = np.mean(temp, axis=-1) #need to add temp alt




###############################################################################

#Import old data
os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL_AM4/reduced_datasets')
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

alt = b['alt'][:]
lat = b['lat'][:]

############################################################################### Get Temperature altitudes
os.chdir('E:/University/University/MSc/Models/Data/CMIP6/gfdl_am4') #Home PC

dataset = Dataset('ta_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
p_temp = np.array(dataset.variables['plev'][:])

#---convert pressure levels to altitude---#

#https://www.mide.com/pages/air-pressure-at-altitude-calculator
#https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
temp_alt_t = np.empty((p_temp.size,1),dtype=float)
temp_alt_p = np.empty((p_temp.size,1),dtype=float)
temp_alt_ts = np.empty((p_temp.size,1),dtype=float)

# Iterate through all of the temp elements (troposphere h < 11km)
i = 0
for item in p_temp:
    newalt = (288.19 - 288.08*((item/101290)**(1/5.256)))/6.49
    temp_alt_t[i] = [newalt]
    i+=1


# Iterate through all of the pressure elements (lower stratosphere 11km < h <25km)
i = 0
for item in p_temp:
    newalt = (1.73 - math.log(item/22650))/0.157
    temp_alt_p[i] = [newalt]
    i+=1


# Iterate through all of the temp elements (upper stratosphere  h > 25km)
i = 0
for item in p_temp:
    newalt = (216.6*((item/2488)**(1/-11.388)) - 141.94)/2.99
    temp_alt_ts[i] = [newalt]
    i+=1


#manually adjust alt and alt_so arrays usinf alt_p and alt_ts
alt_temp = temp_alt_t

###############################################################################

os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL_AM4/reduced_datasets') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GFDL/reduced_datasets') #Uni Laptop

with h5py.File('07.2006_04.2011_GFDL_AM4.h5', 'w') as p:

    p.create_dataset('lat', data=lat)
    p.create_dataset('alt', data=alt)
    p.create_dataset('alt_temp', data=alt_temp) 
    
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
    
    p.create_dataset('lw_alt_lat', data=lw_alt_lat)
    p.create_dataset('iw_alt_lat', data=iw_alt_lat)
    p.create_dataset('cf_alt_lat', data=cf_alt_lat)
    p.create_dataset('temp_alt_lat', data=temp_alt_lat)
    
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















