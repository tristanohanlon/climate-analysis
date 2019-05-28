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
    

#---get comparison altitude from smaller temp and plev arrays---#
i = 0
alt_comp_trop = np.empty(plev.size)

for value in plev:
    altc = (288.19 - 288.08*math.pow((value/101290), (1/5.256))) / 6.49  
    alt_comp_trop[i] = altc
    i+=1
#alt_comp_trop = alt_comp_trop[0:9]  

i = 0
alt_comp_strat = np.empty(plev.size)    

for value in plev:
    altc = ((1.73 - math.log(value/22650)) / 0.157)
    alt_comp_strat[i] = altc
    i+=1
#alt_comp_strat=alt_comp_strat[9:15]  
#alt_comp = np.concatenate((alt_comp_trop, alt_comp_strat), axis = 0)

#---get mass of air path (kg/m^2)---#

#air_density = [] #create empty list
#calculate air density at each pressure layer
#air_density = (plev) / (286.9 * temp)

#air_density = air_density[0:15]

#ap = integrate.trapz(air_density, alt_comp) #mass of air path (kg/m^2)

#---create temperature profile---#
#troposphere from pressure data
#stratosphere assumed from Earth Atmosphere Model 
#https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html

#temp_profile = 288.19 - 6.49*alt

#manually add on 216.69K values to altitudes above 11km to 25km 

###############################################################################
"""
#---create datasets---#

cf = np.vstack((alt, cf)).T
lw = np.vstack((alt, lw)).T
iw = np.vstack((alt, iw)).T
temp = np.vstack((alt, temp_profile)).T
pressure = np.vstack((alt, pressure)).T
    
cf_t = np.vstack((temp_profile, cf[:,1])).T
lw_t = np.vstack((temp_profile, lw[:,1])).T
iw_t = np.vstack((temp_profile, iw[:,1])).T
    
cf_p = np.vstack((pressure[:,1], cf[:,1])).T
lw_p = np.vstack((pressure[:,1], lw[:,1])).T
iw_p = np.vstack((pressure[:,1], iw[:,1])).T


plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()
ax1.plot(cf[:,1],cf[:,0], '-r', label='Fraction Cloud Cover')
ax2.plot(lw[:,1],lw[:,0], '-b', label='Specific Cloud Liquid Water Content')
ax2.plot(iw[:,1],iw[:,0], '--b', label='Specific Cloud Ice Water Content')

#ax.axis('equal')
ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);
ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
          ncol=4, fancybox=True, shadow=True);
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Specific Cloud Liquid and Ice Water Content (kg/kg) x $10^{-4}$')
ax1.set_ylabel('Altitude (km)')
plt.title('Cloud Fraction and Phase vs Altitude GFDL.AM4 2010')
#plt.gca().invert_yaxis()

plt.grid(True)
plt.show()
"""


















