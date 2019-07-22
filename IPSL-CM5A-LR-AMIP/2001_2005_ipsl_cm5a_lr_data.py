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
from scipy import interpolate
"""
#os.chdir('c:/Users/toha006/University/University/MSc/Models/climate-analysis/GFDL-HIRAM-C360/reduced_datasets/backup_reduced_datasets')
os.chdir('//synthesis/e/University/University/MSc/Models/climate-analysis/MRI-CGCM3-AMIP/reduced_datasets/backup_reduced_datasets')
b = h5py.File('2001_2005_mri_cgcm3.h5', 'r')

alt = b['alt'][:]
b.close()
"""

########################################---get variables---########################################

#os.chdir('//synthesis/e/University/University/MSc/Models/Data/CMIP5/mri_cgcm3_amip') #Home PC
os.chdir('D:/MSc/Models/Data/CMIP5/ipsl_cm5a_lr_amip') #ext HDD

f = Dataset('clt_Amon_IPSL-CM5A-LR_amip_r1i1p1_197901-200912.nc', 'r')

#get latitude
lat = np.array(f.variables['lat'][:])

#get total cloud cover keyed to latitude
tcc = np.array(f.variables['clt'][264:324])
tcc = np.mean(tcc, axis = 0)
tcc = np.mean(tcc, axis = -1) / 100
f.close()


#get cloud fraction
f = Dataset('cl_Amon_IPSL-CM5A-LR_amip_r1i1p1_197901-200912.nc', 'r')
cf = np.array(f.variables['cl'][264:324])
cf = np.mean(cf, axis = 0) / 100


#get hybrid pressure levels
plev = np.array(f.variables['lev'][:]) #in hPa
a = np.array(f.variables['ap'][:]) #in hPa
b = np.array(f.variables['b'][:]) #in hPa
#p0 = np.array(f.variables['p0'][:]) #in hPa
f.close()


f = Dataset('ps_Amon_IPSL-CM5A-LR_amip_r1i1p1_197901-200912.nc', 'r')
ps = np.array(f.variables['ps'][264:324]) #in hPa
f.close()

#Convert the hybrid pressure levels to Pa
ps = np.mean(ps, axis = 0)

ps_so = np.mean(ps, axis = -1)
ps_so = np.transpose(ps_so)

ps_so = np.vstack((lat, ps_so)).T #creates a (180,34) array
ps_so = ps_so[ps_so[:,0]>=-70]
ps_so = ps_so[ps_so[:,0]<=-50]
ps_so = ps_so[:,1:] #Split the combined array into just the tccf data, eliminating the first coloumn of latitude
ps_so = np.mean(ps_so, axis = 0)


ps = np.mean(ps, axis = 0)
ps = np.mean(ps, axis = 0)

p = a + b*ps
p = np.array(p / 100) #hPa

p_so = a + b*ps_so
p_so = np.array(p_so / 100) #hPa


#get cloud liquid content
f = Dataset('clw_Amon_IPSL-CM5A-LR_amip_r1i1p1_197901-200912.nc', 'r')
lw = np.array(f.variables['clw'][264:324])
lw = np.mean(lw, axis = 0)
f.close()


#get cloud ice content
f = Dataset('cli_Amon_IPSL-CM5A-LR_amip_r1i1p1_197901-200912.nc', 'r')
iw = np.array(f.variables['cli'][264:324])
iw = np.mean(iw, axis = 0)
f.close()

#get temperature
f = Dataset('ta_Amon_IPSL-CM5A-LR_amip_r1i1p1_197901-200912.nc', 'r')
T = np.array(f.variables['ta'][264:324])
T[T>400] = None
T = np.nanmean(T, axis = 0)
plev_t = np.array(f.variables['plev'][:]) / 100 #in Pa
f.close()

T_so = np.nanmean(T, axis = -1)
T_so = np.transpose(T_so)

T_so = np.hstack((np.vstack(lat), T_so)) #creates a (180,34) array
T_so = T_so[T_so[:,0]>=-70]
T_so = T_so[T_so[:,0]<=-50]
T_so = T_so[:,1:] #Split the combined array into just the tccf data, eliminating the first coloumn of latitude
T_so = np.nanmean(T_so, axis = 0)

T_g = np.nanmean(T, axis = -1)
T_g = np.nanmean(T_g, axis = -1)





###############################################################################
"""
#---convert pressure levels to altitude---#

#https://www.mide.com/pages/air-pressure-at-altitude-calculator
#https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
alt_t = np.empty((p.size,1),dtype=float)
alt_p = np.empty((p.size,1),dtype=float)
alt_ts = np.empty((p.size,1),dtype=float)

# Iterate through all of the temp elements (troposphere h < 11km)
i = 0
for item in p:
    newalt = (288.19 - 288.08*((item/1012.90)**(1/5.256)))/6.49
    alt_t[i] = [newalt]
    i+=1


# Iterate through all of the pressure elements (lower stratosphere 11km < h <25km)
i = 0
for item in p:
    newalt = (1.73 - math.log(item/226.50))/0.157
    alt_p[i] = [newalt]
    i+=1


# Iterate through all of the temp elements (upper stratosphere  h > 25km)
i = 0
for item in p:
    newalt = (216.6*((item/24.88)**(1/-11.388)) - 141.94)/2.99
    alt_ts[i] = [newalt]
    i+=1

sys.exit(0)
#manually adjust alt and alt_so arrays usinf alt_p and alt_ts
alt = alt_t
alt = np.hstack(alt)
"""


#interpolate southern ocean altitudes

f = interpolate.interp1d(p, alt, fill_value="extrapolate", kind = 'cubic')
alt_so = f(p_so)  

f = interpolate.interp1d(p, alt, fill_value="extrapolate", kind = 'cubic')
alt_temp = f(plev_t)

f = interpolate.interp1d(plev_t, T_g,fill_value="extrapolate", kind = 'cubic')
temp_g = f(p)

f = interpolate.interp1d(plev_t, T_so,fill_value="extrapolate", kind = 'cubic')
temp_so = f(p_so)
    
###############################################################################

#---combine arrays---#
# since lw and iw are in kg/kg - need to convert to LWP and IWC in kgm^-2
# Get density levels
# Integrate density with altitude to get air path AP
# multiply lw and iw by AP

air_density = [] #create empty list
#calculate air density at each pressure layer
air_density = ((p * 100) / (286.9 * temp_g))

ap = integrate.trapz(air_density, (alt * 1000)) # gloabl air path

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

temp_alt_lat = np.mean(T, axis = -1) # goes with alt_temp
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
cf_so = np.vstack((alt_so, cf_so)).T


lw_so = np.mean(lw, axis = -1)
lw_so = np.transpose(lw_so)

lw_so = np.hstack((np.vstack(lat), lw_so)) #creates a (180,34) array
lw_so = lw_so[lw_so[:,0]>=-70]
lw_so = lw_so[lw_so[:,0]<=-50]
lw_so = lw_so[:,1:] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
lw_so = np.mean(lw_so, axis = 0)
lw_so = np.vstack((alt_so, lw_so)).T


iw_so = np.mean(iw, axis = -1)
iw_so = np.transpose(iw_so)

iw_so = np.hstack((np.vstack(lat), iw_so)) #creates a (180,34) array
iw_so = iw_so[iw_so[:,0]>=-70]
iw_so = iw_so[iw_so[:,0]<=-50]
iw_so = iw_so[:,1:] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
iw_so = np.mean(iw_so, axis = 0)
iw_so = np.vstack((alt_so, iw_so)).T

#----------------------------#

lwc = lw_so[:,1]

iwc = iw_so[:,1]

lw_frac = (lwc/(lwc+iwc))
iw_frac = (iwc/(lwc+iwc))

lw_frac_so = np.vstack((alt_so, lw_frac * cf_so[:,1])).T
iw_frac_so = np.vstack((alt_so, iw_frac * cf_so[:,1])).T


#----------------------------#
temp_g = np.vstack((alt, temp_g)).T
temp_so = np.vstack((alt_so, temp_so)).T

pressure = np.vstack((alt, p)).T
pressure_so = np.vstack((alt_so, p_so)).T

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




#os.chdir('E:/University/University/MSc/Models/climate-analysis/MRI-CGCM3-AMIP/reduced_datasets') #Home PC
os.chdir('c:/Users/toha006/University/University/MSc/Models/climate-analysis/IPSL-CM5A-LR-AMIP/reduced_datasets') #Home PC
with h5py.File('2001_2005_ipsl_cm5a_lr.h5', 'w') as p:
    
    p.create_dataset('alt', data=alt)
    p.create_dataset('alt_so', data=alt_so)  
    p.create_dataset('alt_temp', data=alt_temp)  
    p.create_dataset('lat', data=lat)  
    p.create_dataset('air_density', data=air_density)  
    p.create_dataset('temp', data=temp_g)
    p.create_dataset('temp_so', data=temp_so)
    p.create_dataset('pressure', data=pressure)
    p.create_dataset('pressure_so', data=pressure_so)
        
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
