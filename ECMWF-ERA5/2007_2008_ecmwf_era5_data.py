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
import sys
import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import h5py
import math
from scipy import integrate
from scipy import interpolate


# Uni Laptop
#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/1979-03.2019_ECMWF_amon_tcc_tciw_tclw.nc', 'r')
# Home PC
dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF-ERA5/1979-03.2019_ECMWF_amon_tcc_tciw_tclw.nc', 'r')

#lon = dataset.variables['longitude'][:] #Extract longitude data
lat = np.array(dataset.variables['latitude'][:]) #Extract latitude data
tcc = np.mean(np.mean(np.array(dataset.variables['tcc'][324:348]), axis = 0), axis = -1) #Extract total cloud cover, keyed to time, lon and lat
tciw = np.mean(np.mean(np.array(dataset.variables['tciw'][324:348]), axis = 0), axis = -1) #Extract ice water path (kg/m^2), keyed to time, lon and lat
tclw = np.mean(np.mean(np.array(dataset.variables['tclw'][324:348]), axis = 0), axis = -1) #Extract liquid water path (kg/m^2), keyed to time, lon and lat

############################################################################### 2007

# Uni Laptop
#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/pressure_levels/2007_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
# Home PC
dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF-ERA5/pressure_levels/2007_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')


p = np.array(dataset.variables['level'][:]) * 100 #Extract pressure (millibars) level data (37 levels)
a_cf = np.mean(np.mean(np.array(dataset.variables['cc'][:]), axis = 0), axis = -1) #Extract fraction of cloud cover, keyed to time, lat and lon (12, 37, 721, 1440)
a_iw = np.mean(np.mean(np.array( dataset.variables['ciwc'][:]), axis = 0), axis = -1) #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
a_lw = np.mean(np.mean(np.array(dataset.variables['clwc'][:]), axis = 0), axis = -1) #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
a_temp = np.mean(np.mean(np.array(dataset.variables['t'][:]), axis = 0), axis = -1) #Extract air temperature profile content (K), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

a_cf_alt_lat = a_cf
a_iw_alt_lat = a_iw
a_lw_alt_lat = a_lw
a_temp_alt_lat = a_temp

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
a_cf_so = np.vstack((lat, a_cf)).T #creates a (721,38) array

#Select latitudes over the southern ocean
a_cf_so = a_cf_so[a_cf_so[:,0]>=-70]
a_cf_so = a_cf_so[a_cf_so[:,0]<=-50]

#Split the combined array into just the cf data, eliminating the first coloumn of latitude
a_cf_so = a_cf_so[:,1:38]

a_cf = np.mean(a_cf, axis=-1) #Average fraction of cloud cover over latitude
a_cf_so = np.mean(a_cf_so, axis=0) #Average southern ocean fraction of cloud cover over latitude


#---------------iw-------------#

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
a_iw_so = np.vstack((lat, a_iw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
a_iw_so = a_iw_so[a_iw_so[:,0]>=-70]
a_iw_so = a_iw_so[a_iw_so[:,0]<=-50]

#Split the combined array into just the iw data, eliminating the first coloumn of latitude
a_iw_so = a_iw_so[:,1:38]

a_iw = np.mean(a_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over latitude
a_iw_so = np.mean(a_iw_so, axis=0) #Average southern ocean specific cloud ice water content (kg/kg) over latitude


#---------------lw-------------#

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
a_lw_so = np.vstack((lat, a_lw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
a_lw_so = a_lw_so[a_lw_so[:,0]>=-70]
a_lw_so = a_lw_so[a_lw_so[:,0]<=-50]

#Split the combined array into just the lw data, eliminating the first coloumn of latitude
a_lw_so = a_lw_so[:,1:38]

a_lw = np.mean(a_lw, axis=-1) #Average specific cloud liquid water content (kg/kg) over latitude
a_lw_so = np.mean(a_lw_so, axis=0) #Average southern ocean specific cloud liquid water content (kg/kg) over latitude


#---------------temperature-------------#


#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
a_temp_so = np.vstack((lat, a_temp)).T #creates a (721,38) array

#Select latitudes over the southern ocean
a_temp_so = a_temp_so[a_temp_so[:,0]>=-70]
a_temp_so = a_temp_so[a_temp_so[:,0]<=-50]

#Split the combined array into just the temp data, eliminating the first coloumn of latitude
a_temp_so = a_temp_so[:,1:38]

a_temp = np.mean(a_temp, axis=-1) #Average air temperature over latitude
a_temp_so = np.mean(a_temp_so, axis=0)  #Average southern ocean air temperature over latitude

############################################################################### 2007

# Uni Laptop
#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/ECMWF/pressure_levels/2008_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
# Home PC
dataset = Dataset('E:/University/University/MSc/Models/Data/ECMWF-ERA5/pressure_levels/2008_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')


b_cf = np.mean(np.mean(np.array(dataset.variables['cc'][:]), axis = 0), axis = -1) #Extract fraction of cloud cover, keyed to time, lat and lon (12, 37, 721, 1440)
b_iw = np.mean(np.mean(np.array( dataset.variables['ciwc'][:]), axis = 0), axis = -1) #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
b_lw = np.mean(np.mean(np.array(dataset.variables['clwc'][:]), axis = 0), axis = -1) #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
b_temp = np.mean(np.mean(np.array(dataset.variables['t'][:]), axis = 0), axis = -1) #Extract air temperature profile content (K), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

b_cf_alt_lat = b_cf
b_iw_alt_lat = b_iw
b_lw_alt_lat = b_lw
b_temp_alt_lat = b_temp

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
b_cf_so = np.vstack((lat, b_cf)).T #creates a (721,38) array

#Select latitudes over the southern ocean
b_cf_so = b_cf_so[b_cf_so[:,0]>=-70]
b_cf_so = b_cf_so[b_cf_so[:,0]<=-50]

#Split the combined array into just the cf data, eliminating the first coloumn of latitude
b_cf_so = b_cf_so[:,1:38]

b_cf = np.mean(b_cf, axis=-1) #Average fraction of cloud cover over latitude
b_cf_so = np.mean(b_cf_so, axis=0) #Average southern ocean fraction of cloud cover over latitude


#---------------iw-------------#

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
b_iw_so = np.vstack((lat, b_iw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
b_iw_so = b_iw_so[b_iw_so[:,0]>=-70]
b_iw_so = b_iw_so[b_iw_so[:,0]<=-50]

#Split the combined array into just the iw data, eliminating the first coloumn of latitude
b_iw_so = b_iw_so[:,1:38]

b_iw = np.mean(b_iw, axis=-1) #Average specific cloud ice water content (kg/kg) over latitude
b_iw_so = np.mean(b_iw_so, axis=0) #Average southern ocean specific cloud ice water content (kg/kg) over latitude


#---------------lw-------------#

#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
b_lw_so = np.vstack((lat, b_lw)).T #creates a (721,38) array

#Select latitudes over the southern ocean
b_lw_so = b_lw_so[b_lw_so[:,0]>=-70]
b_lw_so = b_lw_so[b_lw_so[:,0]<=-50]

#Split the combined array into just the lw data, eliminating the first coloumn of latitude
b_lw_so = b_lw_so[:,1:38]

b_lw = np.mean(b_lw, axis=-1) #Average specific cloud liquid water content (kg/kg) over latitude
b_lw_so = np.mean(b_lw_so, axis=0) #Average southern ocean specific cloud liquid water content (kg/kg) over latitude


#---------------temperature-------------#


#Southern Ocean
#Join the two lists as if they were two columns side by side, into a list of two elements each
b_temp_so = np.vstack((lat, b_temp)).T #creates a (721,38) array

#Select latitudes over the southern ocean
b_temp_so = b_temp_so[b_temp_so[:,0]>=-70]
b_temp_so = b_temp_so[b_temp_so[:,0]<=-50]

#Split the combined array into just the temp data, eliminating the first coloumn of latitude
b_temp_so = b_temp_so[:,1:38]

b_temp = np.mean(b_temp, axis=-1) #Average air temperature over latitude
b_temp_so = np.mean(b_temp_so, axis=0)  #Average southern ocean air temperature over latitude


#################################################################################

#---------------altitude-------------#

alt_t = np.empty((p.size,1),dtype=float)
alt_p = np.empty((p.size,1),dtype=float)
alt_ts = np.empty((p.size,1),dtype=float)

# Iterate through all of the temp elements (troposphere h < 11km)
i = 0
for item in p:
    newalt = (288.19 - 288.08*((item/101290)**(1/5.256)))/6.49
    alt_t[i] = [newalt]
    i+=1


# Iterate through all of the pressure elements (lower stratosphere 11km < h <25km)
i = 0
for item in p:
    newalt = (1.73 - math.log(item/22650))/0.157
    alt_p[i] = [newalt]
    i+=1


# Iterate through all of the temp elements (upper stratosphere  h > 25km)
i = 0
for item in p:
    newalt = (216.6*((item/2488)**(1/-11.388)) - 141.94)/2.99
    alt_ts[i] = [newalt]
    i+=1

sys.exit(0)
#manually adjust alt and alt_so arrays usinf alt_p and alt_ts
alt = alt_t

#################################################################################

#---reduce each dataset---#

############################################################################### cf

combined = a_cf + b_cf

cf = combined / 2

############################################################################### cf_so

combined = a_cf_so + b_cf_so

cf_so = combined / 2

############################################################################### lw

combined = a_lw + b_lw

lw = combined / 2

############################################################################### lw_so

combined = a_lw_so + b_lw_so

lw_so = combined / 2

############################################################################### iw

combined = a_iw + b_iw

iw = combined / 2

############################################################################### iw_so

combined = a_iw_so + b_iw_so

iw_so = combined / 2

############################################################################### temp

combined = a_temp + b_temp

temp_g = combined / 2

############################################################################### temp_so

combined = a_temp_so + b_temp_so

temp_so = combined / 2

############################################################################### cf_alt_lat

combined = a_cf_alt_lat + b_cf_alt_lat

cf_alt_lat = combined / 2

############################################################################### lw_alt_lat

combined = a_lw_alt_lat + b_lw_alt_lat

lw_alt_lat = combined / 2

############################################################################### iw_alt_lat

combined = a_iw_alt_lat + b_iw_alt_lat

iw_alt_lat = combined / 2

############################################################################### temp_alt_lat

combined = a_temp_alt_lat + b_temp_alt_lat

temp_alt_lat = combined / 2

###############################################################################
###############################################################################
###############################################################################
p = p/100

tcc = np.vstack((lat, tcc)).T # Join the two lists as if they were two columns side by side, into a list of two elements each
tclw = np.vstack((lat, tclw)).T
tciw = np.vstack((lat, tciw)).T

#----------------------------#

lw_frac = (tclw[:,1]/(tclw[:,1]+tciw[:,1]))
iw_frac = (tciw[:,1]/(tclw[:,1]+tciw[:,1]))

tclw_frac = np.vstack((lat, lw_frac * tcc[:,1])).T
tciw_frac = np.vstack((lat, iw_frac * tcc[:,1])).T

#----------------------------#

cf_g = np.vstack((alt, cf)).T

lw_g = np.vstack((alt, lw)).T

iw_g = np.vstack((alt, iw)).T

#----------------------------#


lw_frac = (lw_g[:,1]/(lw_g[:,1]+iw_g[:,1]))
iw_frac = (iw_g[:,1]/(lw_g[:,1]+iw_g[:,1]))

lw_frac_g = np.vstack((alt, lw_frac * cf_g[:,1])).T
iw_frac_g = np.vstack((alt, iw_frac * cf_g[:,1])).T


#----------------------------#

cf_so = np.vstack((alt, cf_so)).T

lw_so = np.vstack((alt, lw_so)).T

iw_so = np.vstack((alt, iw_so)).T

#----------------------------#

lw_frac = (lw_so[:,1]/(lw_so[:,1]+iw_so[:,1]))
iw_frac = (iw_so[:,1]/(lw_so[:,1]+iw_so[:,1]))

lw_frac_so = np.vstack((alt, lw_frac * cf_so[:,1])).T
iw_frac_so = np.vstack((alt, iw_frac * cf_so[:,1])).T

#----------------------------#
temp_g = np.vstack((alt, temp_g)).T
temp_so = np.vstack((alt, temp_so)).T

pressure = np.vstack((alt, p)).T

cf_t = np.vstack((temp_g[:,1], cf_g[:,1])).T
cf_t_so = np.vstack((temp_so[:,1], cf_so[:,1])).T
lw_t = np.vstack((temp_g[:,1], lw_g[:,1])).T
lw_t_so = np.vstack((temp_so[:,1], lw_so[:,1])).T
iw_t = np.vstack((temp_g[:,1], iw_g[:,1])).T
iw_t_so = np.vstack((temp_so[:,1], iw_so[:,1])).T

#----------------------------#

lw_frac_t = np.vstack((temp_g[:,1], lw_frac_g[:,1])).T
lw_frac_t_so = np.vstack((temp_so[:,1], lw_frac_so[:,1])).T
iw_frac_t = np.vstack((temp_g[:,1], iw_frac_g[:,1])).T
iw_frac_t_so = np.vstack((temp_so[:,1], iw_frac_so[:,1])).T


fig, ax1 = plt.subplots()
ax1.plot(lw_frac_t_so[:,0],lw_frac_t_so[:,1], '-r', label='SO')
ax1.plot(lw_frac_t[:,0],lw_frac_t[:,1], '-b', label='G')
#----------------------------#




os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF-ERA5/reduced_datasets')
with h5py.File('2007_2008_ecmwf_era5.h5', 'w') as p:
    
    p.create_dataset('alt', data=alt)
    p.create_dataset('lat', data=lat)  
    p.create_dataset('temp', data=temp_g)
    p.create_dataset('temp_so', data=temp_so)
    p.create_dataset('pressure', data=pressure)
        
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('tclw', data=tclw)
    p.create_dataset('tciw', data=tciw)
    p.create_dataset('tclw_frac', data=tclw_frac)
    p.create_dataset('tciw_frac', data=tciw_frac)
    
    p.create_dataset('cf', data=cf_g)
    p.create_dataset('lw', data=lw_g)
    p.create_dataset('iw', data=iw_g)
    p.create_dataset('lw_frac', data=lw_frac_g)
    p.create_dataset('iw_frac', data=iw_frac_g)
    
    p.create_dataset('cf_so', data=cf_so)
    p.create_dataset('lw_so', data=lw_so)
    p.create_dataset('iw_so', data=iw_so)
    p.create_dataset('lw_frac_so', data=lw_frac_so)
    p.create_dataset('iw_frac_so', data=iw_frac_so)  
    
    p.create_dataset('temp_alt_lat', data=temp_alt_lat)  
    p.create_dataset('cf_alt_lat', data=cf_alt_lat)  
    p.create_dataset('lw_alt_lat', data=lw_alt_lat)
    p.create_dataset('iw_alt_lat', data=iw_alt_lat)
    
    p.create_dataset('cf_t', data=cf_t)
    p.create_dataset('cf_t_so', data=cf_t_so)

    p.create_dataset('lw_t', data=lw_t)
    p.create_dataset('lw_t_so', data=lw_t_so)
    p.create_dataset('lw_frac_t', data=lw_frac_t)
    p.create_dataset('lw_frac_t_so', data=lw_frac_t_so)

    p.create_dataset('iw_t', data=iw_t)
    p.create_dataset('iw_t_so', data=iw_t_so)
    p.create_dataset('iw_frac_t', data=iw_frac_t)
    p.create_dataset('iw_frac_t_so', data=iw_frac_t_so)
 
 
    p.close()
