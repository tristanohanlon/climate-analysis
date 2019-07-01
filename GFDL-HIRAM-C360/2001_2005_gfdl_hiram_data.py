# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 18:40:37 2019

@author: Tristan O'Hanlon

This will create a dataset of time averaged global total cloud cover with latitude.
Time period is from 01.1980 - 12.2014

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


os.chdir('//synthesis/e/University/University/MSc/Models/Data/CMIP5/gfdl_hiram_c360') #Home PC
#os.chdir('D:/MSc/Models/Data/CMIP6/cesm2.1_cam6') #ext HDD

f = Dataset('clt_Amon_GFDL-HIRAM-C360_amip_r1i1p1_199901-200312.nc', 'r')

#get latitude
lat = np.array(f.variables['lat'][:])

#get total cloud cover keyed to latitude
tcc1 = np.array(f.variables['clt'][24:])

f = Dataset('clt_Amon_GFDL-HIRAM-C360_amip_r1i1p1_200401-200812.nc', 'r')
tcc2 = np.array(f.variables['clt'][:24])
tcc = np.concatenate((tcc1, tcc2))

tcc = np.mean(tcc, axis = 0)
tcc = np.mean(tcc, axis = -1) / 100


#get cloud fraction
f = Dataset('cl_Amon_GFDL-HIRAM-C360_amip_r1i1p1_199901-200312.nc', 'r')
cf1 = np.array(f.variables['cl'][24:])
f = Dataset('cl_Amon_GFDL-HIRAM-C360_amip_r1i1p1_200401-200812.nc', 'r')
cf2 = np.array(f.variables['cl'][:24])
cf = np.concatenate((cf1, cf2))

cf = np.mean(cf, axis = 0)


#get hybrid pressure levels
plev = np.array(f.variables['lev'][:]) #in hPa
a = np.array(f.variables['a'][:]) #in hPa
b = np.array(f.variables['b'][:]) #in hPa
p0 = np.array(f.variables['p0'][:]) #in hPa
f = Dataset('ps_Amon_GFDL-HIRAM-C360_amip_r1i1p1_199901-200312.nc', 'r')
ps1 = np.array(f.variables['ps'][24:]) #in hPa
f = Dataset('ps_Amon_GFDL-HIRAM-C360_amip_r1i1p1_200401-200812.nc', 'r')
ps2 = np.array(f.variables['ps'][:24])
ps = np.concatenate((ps1, ps2))

#Convert the hybrid pressure levels to Pa
ps = np.mean(ps, axis = 0)
ps = np.mean(ps, axis = 0)
ps = np.mean(ps, axis = 0)

p = a*p0 + b*ps
p = np.array(p)


#get cloud liquid content
f = Dataset('clw_Amon_GFDL-HIRAM-C360_amip_r1i1p1_199901-200312.nc', 'r')
lw1 = np.array(f.variables['clw'][24:])
f = Dataset('clw_Amon_GFDL-HIRAM-C360_amip_r1i1p1_200401-200812.nc', 'r')
lw2 = np.array(f.variables['clw'][:24])
lw = np.concatenate((lw1, lw2))
lw = np.mean(lw, axis = 0)


#get cloud ice content
f = Dataset('cli_Amon_GFDL-HIRAM-C360_amip_r1i1p1_199901-200312.nc', 'r')
iw1 = np.array(f.variables['cli'][24:])
f = Dataset('cli_Amon_GFDL-HIRAM-C360_amip_r1i1p1_200401-200812.nc', 'r')
iw2 = np.array(f.variables['cli'][:24])
iw = np.concatenate((iw1, iw2))
iw = np.mean(iw, axis = 0)

#get temperature
f = Dataset('ta_Amon_GFDL-HIRAM-C360_amip_r1i1p1_199901-200312.nc', 'r')
T1 = np.array(f.variables['ta'][24:])
f = Dataset('ta_Amon_GFDL-HIRAM-C360_amip_r1i1p1_200401-200812.nc', 'r')
T2 = np.array(f.variables['ta'][:24])
T = np.concatenate((T1, T2))
T = np.mean(T, axis = 0)
T[T>400] = None
plev = np.array(f.variables['plev'][:]) #in Pa



###############################################################################

#---convert pressure levels to altitude---#

#https://www.mide.com/pages/air-pressure-at-altitude-calculator
#https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
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



###############################################################################

alt_t = np.empty((plev.size,1),dtype=float)
alt_p = np.empty((plev.size,1),dtype=float)
alt_ts = np.empty((plev.size,1),dtype=float)

# Iterate through all of the temp elements (troposphere h < 11km)
i = 0
for item in plev:
    newalt = (288.19 - 288.08*((item/101290)**(1/5.256)))/6.49
    alt_t[i] = [newalt]
    i+=1


# Iterate through all of the pressure elements (lower stratosphere 11km < h <25km)
i = 0
for item in plev:
    newalt = (1.73 - math.log(item/22650))/0.157
    alt_p[i] = [newalt]
    i+=1


# Iterate through all of the temp elements (upper stratosphere  h > 25km)
i = 0
for item in plev:
    newalt = (216.6*((item/2488)**(1/-11.388)) - 141.94)/2.99
    alt_ts[i] = [newalt]
    i+=1

temp_g = np.mean(T, axis = -1)
temp_g = np.mean(temp_g, axis = -1)
temp_g = np.vstack(temp_g)
temp_g = np.hstack((np.vstack(alt_t), temp_g))

#-----------------------#
alt_temp = 288.14 - 6.49 * alt    
alt_ts = 141.89 + 2.99 * alt
    
    
###############################################################################

#---combine arrays---#
# since lw and iw are in kg/kg - need to convert to LWP and IWC in kgm^-2
# Get density levels
# Integrate density with altitude to get air path AP
# multiply lw and iw by AP

alt_t = np.hstack(alt*1000)
pressure = np.hstack((np.vstack(alt), np.vstack(p)))

air_density = [] #create empty list
#calculate air density at each pressure layer
air_density = (pressure[:,1]) / (286.9 * alt_temp[:,0])

ap = integrate.trapz(air_density, alt_t)

tclw = np.mean(lw , axis = 0)
tclw = np.mean(tclw  , axis = -1) * ap

tciw = np.mean(iw , axis = 0)
tciw = np.mean(tciw , axis = -1) * ap

tcc = np.vstack((lat, tcc)).T # Join the two lists as if they were two columns side by side, into a list of two elements each
tclw = np.vstack((lat, tclw)).T
tciw = np.vstack((lat, tciw)).T

#----------------------------#

lwc = np.mean(lw , axis = 0)
lwc = np.mean(lwc , axis = -1)

iwc = np.mean(iw , axis = 0)
iwc = np.mean(iwc , axis = -1)

lw_frac = (lwc/(lwc+iwc))
iw_frac = (iwc/(lwc+iwc))

tclw_frac = np.vstack((lat, lw_frac * tcc[:,1])).T
tciw_frac = np.vstack((lat, iw_frac * tcc[:,1])).T

#----------------------------#
alt = np.hstack(alt)

cf_g = np.mean(cf, axis = -1)
cf_g = np.mean(cf_g, axis = -1)
cf_g = np.vstack((alt, cf_g)).T

lw_g = np.mean(lw, axis = -1)
lw_g = np.mean(lw_g, axis = -1)
lw_g = np.vstack((alt, lw_g)).T

iw_g = np.mean(iw, axis = -1)
iw_g = np.mean(iw_g, axis = -1)
iw_g = np.vstack((alt, iw_g)).T

#----------------------------#

lwc = np.mean(lw , axis = -1)
lwc = np.mean(lwc , axis = -1)

iwc = np.mean(iw , axis = -1)
iwc = np.mean(iwc , axis = -1)

lw_frac = (lwc/(lwc+iwc))
iw_frac = (iwc/(lwc+iwc))

lw_frac_g = np.vstack((alt, lw_frac * cf_g[:,1])).T
iw_frac_g = np.vstack((alt, iw_frac * cf_g[:,1])).T

#----------------------------#

temp_alt_lat = np.mean(T, axis = -1)
temp_alt_lat[temp_alt_lat>400] = None
cf_alt_lat = np.mean(cf, axis = -1)
lw_alt_lat = np.mean(lw, axis = -1)
iw_alt_lat = np.mean(iw, axis = -1)

#Select Southern ocean Latitudes
cf_so = np.mean(cf, axis = -1)
cf_so = np.transpose(cf_so)

cf_so = np.hstack((np.vstack(lat), cf_so)) #creates a (180,34) array
cf_so = cf_so[cf_so[:,0]>=-70]
cf_so = cf_so[cf_so[:,0]<=-50]
cf_so = cf_so[:,1:] #Split the combined array into just the tccf data, eliminating the first coloumn of latitude
cf_so = np.mean(cf_so, axis = 0)
cf_so = np.vstack((alt, cf_so)).T


lw_so = np.mean(lw, axis = -1)
lw_so = np.transpose(lw_so)

lw_so = np.hstack((np.vstack(lat), lw_so)) #creates a (180,34) array
lw_so = lw_so[lw_so[:,0]>=-70]
lw_so = lw_so[lw_so[:,0]<=-50]
lw_so = lw_so[:,1:] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
lw_so = np.mean(lw_so, axis = 0)
lw_so = np.vstack((alt, lw_so)).T


iw_so = np.mean(iw, axis = -1)
iw_so = np.transpose(iw_so)

iw_so = np.hstack((np.vstack(lat), iw_so)) #creates a (180,34) array
iw_so = iw_so[iw_so[:,0]>=-70]
iw_so = iw_so[iw_so[:,0]<=-50]
iw_so = iw_so[:,1:] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
iw_so = np.mean(iw_so, axis = 0)
iw_so = np.vstack((alt, iw_so)).T

#----------------------------#

lwc = np.mean(lw_so , axis = -1)
lwc = np.mean(lwc , axis = -1)

iwc = np.mean(iw_so , axis = -1)
iwc = np.mean(iwc , axis = -1)

lw_frac = (lwc/(lwc+iwc))
iw_frac = (iwc/(lwc+iwc))

lw_frac_so = np.vstack((alt, lw_frac * cf_so[:,1])).T
iw_frac_so = np.vstack((alt, iw_frac * cf_so[:,1])).T


#----------------------------#
temp_g = np.hstack((np.vstack(alt), alt_temp))


cf_t = np.vstack((temp_g[:,1], cf_g[:,1])).T
cf_t_so = np.vstack((temp_g[:,1], cf_so[:,1])).T
lw_t = np.vstack((temp_g[:,1], lw_g[:,1])).T
lw_t_so = np.vstack((temp_g[:,1], lw_so[:,1])).T
iw_t = np.vstack((temp_g[:,1], iw_g[:,1])).T
iw_t_so = np.vstack((temp_g[:,1], iw_so[:,1])).T

#----------------------------#

lw_frac_t = np.vstack((temp_g[:,1], lw_frac_g[:,1])).T
lw_frac_t_so = np.vstack((temp_g[:,1], lw_frac_so[:,1])).T
iw_frac_t = np.vstack((temp_g[:,1], iw_frac_g[:,1])).T
iw_frac_t_so = np.vstack((temp_g[:,1], iw_frac_so[:,1])).T

#----------------------------#




os.chdir('c:/Users/tristan/University/University/MSc/Models/climate-analysis/GFDL-HIRAM-C360/reduced_datasets') #Home PC
with h5py.File('2001_2005_gfdl_hiram.h5', 'w') as p:
    
    p.create_dataset('alt', data=alt)
    p.create_dataset('alt_temp', data=alt_t)  
    p.create_dataset('lat', data=lat)  
    p.create_dataset('air_density', data=air_density)  
    p.create_dataset('temp', data=temp_g)
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





"""
fig, ax1 = plt.subplots()
ax1.plot(lat,tcc, '-r', label='Total Cloud Fraction')
"""