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

"""
#---get cf---#

#os.chdir('E:/University/University/MSc/Models/Data/CMIP6/mri_esm2/amip/cf') #Home PC
os.chdir('//Synthesis/E/University/University/MSc/Models/Data/CMIP6/mri_esm2/cf') #Home PC

f = Dataset('cl_Amon_MRI-ESM2-0_amip_r1i1p1f1_gn_199901-200812.nc', 'r')
cf1 = np.array(f.variables['cl'][90:])
f = Dataset('cl_Amon_MRI-ESM2-0_amip_r1i1p1f1_gn_200901-201412.nc', 'r')
cf2 = np.array(f.variables['cl'][:24])
cf3 = np.array(f.variables['cl'][25:28])

cf = np.concatenate((cf1, cf2, cf3))
cf = np.mean(cf, axis = 0) / 100 # average over time

#---get lw---#

#os.chdir('E:/University/University/MSc/Models/Data/CMIP6/mri_esm2/amip/tclw') #Home PC
os.chdir('//synthesis/e/University/University/MSc/Models/Data/CMIP6/mri_esm2/tclw') #Home PC

f = Dataset('clw_Amon_MRI-ESM2-0_amip_r1i1p1f1_gn_199901-200812.nc', 'r')
lw1 = np.array(f.variables['clw'][90:])
f = Dataset('clw_Amon_MRI-ESM2-0_amip_r1i1p1f1_gn_200901-201412.nc', 'r')
lw2 = np.array(f.variables['clw'][:24])
f = Dataset('clw_Amon_MRI-ESM2-0_amip_r1i1p1f1_gn_200901-201412.nc', 'r')
lw3 = np.array(f.variables['clw'][25:28])
    
lw = np.concatenate((lw1, lw2, lw3))
lw = np.mean(lw, axis = 0)

#---get iw---#

#os.chdir('E:/University/University/MSc/Models/Data/CMIP6/mri_esm2/amip/tciw') #Home PC
os.chdir('//synthesis/e/University/University/MSc/Models/Data/CMIP6/mri_esm2/tciw') #Home PC

f = Dataset('cli_Amon_MRI-ESM2-0_amip_r1i1p1f1_gn_199901-200812.nc', 'r')
iw1 = np.array(f.variables['cli'][90:])
f = Dataset('cli_Amon_MRI-ESM2-0_amip_r1i1p1f1_gn_200901-201412.nc', 'r')
iw2 = np.array(f.variables['cli'][:24])
f = Dataset('cli_Amon_MRI-ESM2-0_amip_r1i1p1f1_gn_200901-201412.nc', 'r')
iw3 = np.array(f.variables['cli'][25:28])
 
iw = np.concatenate((iw1, iw2, iw3))
iw = np.mean(iw, axis = 0)

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

#Average over altitude

lw = np.mean(lw, axis=0)
iw = np.mean(iw, axis=0)

#Average over longitude

lw = np.mean(lw, axis=-1)
iw = np.mean(iw, axis=-1)
"""


###############################################################################

#Import old data
#os.chdir('//synthesis/e/University/University/MSc/Models/climate-analysis/MRI-ESM2-AMIP/reduced_datasets')
os.chdir('c:/Users/tristan/University/University/MSc/Models/climate-analysis/MRI-ESM2-AMIP/reduced_datasets')
b = h5py.File('07.2006_04.2011_mri_esm2.h5', 'r')

tcc = b['tcc'][:] # 0-1
tclw = b['tclw'][:] / 100
tciw = b['tciw'][:] / 100
lat = b['lat'][:]

tclw_frac = b['tclw_frac'][:]
tciw_frac = b['tciw_frac'][:]

cf_alt_lat = b['cf_alt_lat'][:]
lw_alt_lat = b['lw_alt_lat'][:] #kg/kg
iw_alt_lat = b['iw_alt_lat'][:] #kg/kg
temp_alt_lat = b['temp_alt_lat'][:] #kg/kg
alt = b['alt'][:] 
 
cf = b['cf'][:] # 0-1
lw = b['lw'][:] #kg/kg
iw = b['iw'][:] #kg/kg
temp = b['temp'][:] #K
pressure = b['pressure'][:] #hPa

cf_t = b['cf_t'][:] # 0-1
lw_t = b['lw_t'][:] #kg/kg
iw_t = b['iw_t'][:] #kg/kg


cf_so = b['cf_so'][:] # 0-1
lw_so = b['lw_so'][:] #kg/kg
iw_so = b['iw_so'][:] #kg/kg

cf_t = b['cf_t'][:] # 0-1
lw_t = b['lw_t'][:] #kg/kg
iw_t = b['iw_t'][:] #kg/kg


cf_t_so = b['cf_t_so'][:] # 0-1
lw_t_so = b['lw_t_so'][:] #kg/kg
iw_t_so = b['iw_t_so'][:] #kg/kg

lw_frac_so = b['lw_frac_so'][:]
iw_frac_so = b['iw_frac_so'][:]

#--------------------------#

lw_frac_t = b['lw_frac_t'][:]
lw_frac_t_so = b['lw_frac_t_so'][:]
iw_frac_t = b['iw_frac_t'][:]
iw_frac_t_so = b['iw_frac_t_so'][:]

lw_frac_g = b['lw_frac'][:]
iw_frac_g = b['iw_frac'][:]

b.close()
#--------------------------#
"""
lw_frac = (lw[:,1]/(lw[:,1]+iw[:,1]))
iw_frac = (iw[:,1]/(lw[:,1]+iw[:,1]))

lw_frac_g = np.vstack((alt, lw_frac * cf[:,1])).T
iw_frac_g = np.vstack((alt, iw_frac * cf[:,1])).T

#--------------------------#

lw_frac = (lw_so[:,1]/(lw_so[:,1]+iw_so[:,1]))
iw_frac = (iw_so[:,1]/(lw_so[:,1]+iw_so[:,1]))

lw_frac_so = np.vstack((alt, lw_frac * cf_so[:,1])).T
iw_frac_so = np.vstack((alt, iw_frac * cf_so[:,1])).T

#--------------------------#

lw_frac_t = np.vstack((temp[:,1], lw_frac_g[:,1])).T
lw_frac_t_so = np.vstack((temp[:,1], lw_frac_so[:,1])).T
iw_frac_t = np.vstack((temp[:,1], iw_frac_g[:,1])).T
iw_frac_t_so = np.vstack((temp[:,1], iw_frac_so[:,1])).T
"""
#--------------------------#


fig, ax = plt.subplots()
ax.plot(tclw_frac[:,0], tclw_frac[:,1], '-b')
ax.plot(tciw_frac[:,0], tciw_frac[:,1], '--b')



###############################################################################

os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/MRI-ESM2-AMIP/reduced_datasets') # Home PC

with h5py.File('07.2006_04.2011_mri_esm2.h5', 'w') as p:

    p.create_dataset('lat', data=lat)
    p.create_dataset('alt', data=alt)
    p.create_dataset('temp', data=temp)
    p.create_dataset('pressure', data=pressure)
    
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('tclw', data=tclw)
    p.create_dataset('tciw', data=tciw)  
    p.create_dataset('tclw_frac', data=tclw_frac)
    p.create_dataset('tciw_frac', data=tciw_frac)  
    
    p.create_dataset('cf', data=cf)
    p.create_dataset('lw', data=lw)
    p.create_dataset('iw', data=iw)
    p.create_dataset('lw_frac', data=lw_frac_g)
    p.create_dataset('iw_frac', data=iw_frac_g)
    
    p.create_dataset('cf_so', data=cf_so)
    p.create_dataset('lw_so', data=lw_so)
    p.create_dataset('iw_so', data=iw_so)
    p.create_dataset('lw_frac_so', data=lw_frac_so)
    p.create_dataset('iw_frac_so', data=iw_frac_so)  
    
    p.create_dataset('lw_alt_lat', data=lw_alt_lat)
    p.create_dataset('iw_alt_lat', data=iw_alt_lat)
    p.create_dataset('cf_alt_lat', data=cf_alt_lat)
    p.create_dataset('temp_alt_lat', data=temp_alt_lat)
    
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












