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


#---get latitude and cf---#
#os.chdir('E:/University/University/MSc/Models/Data/CMIP6/gfdl_am4') #Home PC
os.chdir('//Synthesis/E/University/University/MSc/Models/Data/CMIP6/cesm2.1_cam6/') #Home PC

f = Dataset('clt_Amon_CESM2_amip_r1i1p1f1_gn_195001-201412.nc', 'r') # 780 months
tcc = np.array(f.variables['CLDTOT'][612:672])

lat = f.variables['lat'][:]
lat = np.array(lat)

tcc = np.mean(tcc, axis = 0) # average over time
tcc = np.mean(tcc, axis = -1) # average over longitude
f.close()

#---get cf---#

#os.chdir('E:/University/University/MSc/Models/Data/CMIP6/mri_esm2/amip/cf') #Home PC

f = Dataset('cl_Amon_CESM2_amip_r1i1p1f1_gn_195001-201412.nc', 'r') # 780 months
cf = np.array(f.variables['CLOUD'][612:672])
cf = np.mean(cf, axis = 0) # average over time

plev = np.array(f.variables['lev'][:]) #in hPa
a = np.array(f.variables['hyam'][:]) #in hPa
b = np.array(f.variables['hybm'][:]) #in hPa
p0 = np.array(f.variables['P0'][:]) #in hPa
f.close()



#---get surface pressure---#

f = Dataset('ps_Amon_CESM2_amip_r1i1p1f1_gn_195001-201412.nc', 'r') # 780 months
ps = np.array(f.variables['PS'][612:672])

f.close()

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

p = a*p0 + b*ps
p = np.array(p / 100) #hPa

p_so = a*p0 + b*ps_so
p_so = np.array(p_so / 100) #hPa



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

"""
os.chdir('E:/University/University/MSc/Models/climate-analysis/CESM2-CAM6-AMIP/reduced_datasets/backup_reduced_datasets')
b = h5py.File('2001_2005_cesm2_cam6.h5', 'r')

alt = b['alt'][:]
b.close()
"""

#interpolate southern ocean altitudes

f = interpolate.interp1d(p, alt, fill_value="extrapolate", kind = 'cubic')
alt_so = f(p_so) 

###############################################################################


#---get lw---#

#os.chdir('E:/University/University/MSc/Models/Data/CMIP6/mri_esm2/amip/tclw') #Home PC
os.chdir('//synthesis/e/University/University/MSc/Models/Data/CMIP6/cesm2.1_cam6') #Home PC

f = Dataset('clw_Amon_CESM2_amip_r1i1p1f1_gn_195001-201412.nc', 'r') # 780 months
lw = np.array(f.variables['CLDLIQ'][612:672])
lw = np.mean(lw, axis = 0)
f.close()

#---get iw---#

#os.chdir('E:/University/University/MSc/Models/Data/CMIP6/mri_esm2/amip/tciw') #Home PC
os.chdir('//synthesis/e/University/University/MSc/Models/Data/CMIP6/cesm2.1_cam6') #Home PC

f = Dataset('cli_Amon_CESM2_amip_r1i1p1f1_gn_195001-201412.nc', 'r') # 780 months
iw = np.array(f.variables['CLDICE'][612:672])
iw = np.mean(iw, axis = 0)
f.close()


#---get temp---#

#os.chdir('E:/University/University/MSc/Models/Data/CMIP6/mri_esm2/amip/') #Home PC
os.chdir('//synthesis/e/University/University/MSc/Models/Data/CMIP6/cesm2.1_cam6') #Home PC

f = Dataset('ta_Amon_CESM2_amip_r1i1p1f1_gn_195001-201412.nc', 'r') # 780 months
T = np.array(f.variables['T'][612:672])
T[T>400] = None
T = np.nanmean(T, axis = 0)
f.close()

T_so = np.nanmean(T, axis = -1)
T_so = np.transpose(T_so)

T_so = np.hstack((np.vstack(lat), T_so)) #creates a (180,34) array
T_so = T_so[T_so[:,0]>=-70]
T_so = T_so[T_so[:,0]<=-50]
T_so = T_so[:,1:] #Split the combined array into just the tccf data, eliminating the first coloumn of latitude
temp_so = np.nanmean(T_so, axis = 0)

T_g = np.nanmean(T, axis = -1)
temp_g = np.nanmean(T_g, axis = -1)
    
    
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
tclw = np.mean(tclw  , axis = -1) * -ap

tciw = np.mean(iw , axis = 0)
tciw = np.mean(tciw , axis = -1) * -ap

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

lwc = np.mean(lw_so , axis = -1)
lwc = np.mean(lwc , axis = -1)

iwc = np.mean(iw_so , axis = -1)
iwc = np.mean(iwc , axis = -1)

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



os.chdir('E:/University/University/MSc/Models/climate-analysis/CESM2-CAM6-AMIP/reduced_datasets') #Home PC
with h5py.File('2001_2005_cesm2_cam6.h5', 'w') as p:
    
    p.create_dataset('alt', data=alt)
    p.create_dataset('alt_so', data=alt_so)  
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