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
os.chdir('//synthesis/e/University/University/MSc/Models/Data/CMIP6/gfdl_am4') #Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/Data/GFDL/') #Uni Laptop
dataset = Dataset('cl_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
cf = dataset.variables['cl'][:] #Extract Cloud Fraction, keyed to time, level, lon and lat (420, 33, 180, 288)

dataset = Dataset('clw_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
lw = dataset.variables['clw'][:] #Extract Mass Fraction of Cloud Liquid Water (kg/kg), keyed to time, level, lon and lat (420, 33, 180, 288)

dataset = Dataset('cli_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
iw = dataset.variables['cli'][:] #Extract Mass Fraction of Cloud Ice (kg/kg), keyed to time, level, lon and lat (420, 33, 180, 288)
end = time.time()
print('Importing data took:', end - start, 's')

#----------------------------#

#Select the months from 07.2006 to 04.2011
cf = np.take(cf, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 
lw = np.take(lw, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 
iw = np.take(iw, [319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376], axis=0) 

#Average over time
cf = np.mean(cf, axis=1)
lw = np.mean(lw, axis=1)
iw = np.mean(iw, axis=1)

###############################################################################

#Import old data
os.chdir('//synthesis/e/University/University/MSc/Models/climate-analysis/GFDL-AM4-AMIP/reduced_datasets')
b = h5py.File('07.2006_04.2011_GFDL_AM4.h5', 'r')

tcc = b['tcc'][:] # 0-1
tclw = b['tclw'][:] #kgm^-2
tciw = b['tciw'][:] #kgm^-2
cf_alt_lat = b['cf_alt_lat'][:]
lw_alt_lat = b['lw_alt_lat'][:] #kg/kg
iw_alt_lat = b['iw_alt_lat'][:] #kg/kg
temp_alt_lat = b['temp_alt_lat'][:] #kg/kg
alt = b['alt'][:] 
alt_temp = b['alt_temp'][:] 
 
lat = b['lat'][:]


cf = b['cf'][:] # 0-1
lw = b['lw'][:] #kg/kg
iw = b['iw'][:] #kg/kg
temp = b['temp'][:] #K
pressure = b['pressure'][:] #hPa

gfdl4_tcc_temp_g = b['cf_t'][:] # 0-1
gfdl4_tclw_temp_g = b['lw_t'][:] #kg/kg
gfdl4_tciw_temp_g = b['iw_t'][:] #kg/kg


cf_so = b['cf_so'][:23] # 0-1
lw_so = b['lw_so'][:19] #kg/kg
iw_so = b['iw_so'][:23] #kg/kg
temp_so = b['temp_so'][:] #K
pressure_so = b['pressure_so'][:] #hPa

cf_t = b['cf_t'][:] # 0-1
lw_t = b['lw_t'][:] #kg/kg
iw_t = b['iw_t'][:] #kg/kg


cf_t_so = b['cf_t_so'][:] # 0-1
lw_t_so = b['lw_t_so'][:] #kg/kg
iw_t_so = b['iw_t_so'][:] #kg/kg


lw_frac = (lw/(lw+iw))
iw_frac = (iw/(lw+iw))

tclw_frac = lw_frac * tcc[:,1]
tciw_frac = iw_frac * tcc[:,1]


fig, ax = plt.subplots()
ax.plot(lat, tclw_frac_tcc, '-k')
ax.plot(lat, tclw_frac, '-b')
ax.plot(lat, tciw_frac_tcc, '--k')
ax.plot(lat, tciw_frac, '--b')



###############################################################################

os.chdir('//synthesis/e/University/University/MSc/Models/climate-analysis/GFDL-AM4-AMIP/reduced_datasets') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GFDL/reduced_datasets') #Uni Laptop

with h5py.File('07.2006_04.2011_GFDL_AM41.h5', 'w') as p:

    p.create_dataset('lat', data=lat)
    p.create_dataset('alt', data=alt)
    p.create_dataset('alt_temp', data=alt_temp) 
    
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('tclw', data=tclw)
    p.create_dataset('tciw', data=tciw)  
    p.create_dataset('tclw_frac', data=tclw_frac)
    p.create_dataset('tciw_frac', data=tciw_frac)  
    
    p.create_dataset('cf', data=cf)
    p.create_dataset('cf_so', data=cf_so)
    p.create_dataset('lw', data=lw)
    p.create_dataset('lw_so', data=lw_so)
    p.create_dataset('iw', data=iw)
    p.create_dataset('iw_so', data=iw_so)

    p.create_dataset('temp', data=temp)
    p.create_dataset('temp_so', data=temp_so)
    p.create_dataset('pressure', data=pressure)
    p.create_dataset('pressure_so', data=pressure_so)
    
    p.create_dataset('lw_alt_lat', data=lw_alt_lat)
    p.create_dataset('iw_alt_lat', data=iw_alt_lat)
    p.create_dataset('cf_alt_lat', data=cf_alt_lat)
    p.create_dataset('temp_alt_lat', data=temp_alt_lat)
    
    p.create_dataset('cf_t', data=cf_t)
    p.create_dataset('cf_t_so', data=cf_t_so)
    p.create_dataset('lw_t', data=lw_t)
    p.create_dataset('lw_t_so', data=lw_t_so)
    p.create_dataset('iw_t', data=iw_t)
    p.create_dataset('iw_t_so', data=iw_t_so)

    p.close()












