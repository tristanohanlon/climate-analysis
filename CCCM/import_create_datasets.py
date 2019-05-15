# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon
"""

#---Import Latitude Variables---#

import h5py
import os
import numpy as np
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM') #Uni Laptop
#os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CCCM') #Laptop


g = h5py.File('raw_datasets/2011_CCCM_tcc_lat.h5', 'r')
tcc = g['tcc'][:]

f = h5py.File('raw_datasets/2011_CCCM_tclw_lat.h5', 'r')
tclw = f['tclw'][:]

h = h5py.File('raw_datasets/2011_CCCM_tciw_lat.h5', 'r')
tciw = h['tciw'][:]


b = h5py.File('raw_datasets/2011_CCCM_profile_variables.h5', 'r')
lat = b['lat'][:]
alt = b['alt'][:]
alt_c = b['alt_c'][:]
alt_t = b['alt_t'][:]
air_density_g = b['air_density_g'][:]
air_density_so = b['air_density_so'][:]
temp = b['temp_g_alt'][:]
temp_so = b['temp_so_alt'][:]
pressure = b['pressure_g_alt'][:]
pressure_so = b['pressure_so_alt'][:]

c = h5py.File('raw_datasets/2011_CCCM_tcc_alt.h5', 'r')
cf = c['cf'][:]
cf_so = c['cf_so'][:]

d = h5py.File('raw_datasets/2011_CCCM_tclw_alt.h5', 'r')
lw = d['lw'][:]
lw_so = d['lw_so'][:]

e = h5py.File('raw_datasets/2011_CCCM_tciw_alt.h5', 'r')
iw = e['iw'][:]
iw_so = e['iw_so'][:]

#---Fit temperture data to cloud fraction and phase---#

cf_scale = cf[12:109]
cf_scale_so = cf_so[12:109]

t_scale_c = temp[:,1]
t_scale_c = t_scale_c[36:133]

t_scale_liw = temp[:,1]
t_scale_liw = t_scale_liw[0:137]

cf_t = np.vstack((t_scale_c, cf_scale[:,1])).T
cf_t_so = np.vstack((t_scale_c, cf_scale_so[:,1])).T
lw_t = np.vstack((t_scale_liw, lw[:,1])).T
lw_t_so = np.vstack((t_scale_liw, lw_so[:,1])).T
iw_t = np.vstack((t_scale_liw, iw[:,1])).T
iw_t_so = np.vstack((t_scale_liw, iw_so[:,1])).T

###############################################################################
    

import h5py
import os
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Uni Laptop

with h5py.File('2011_CCCM.h5', 'w') as p:
    p.create_dataset('lat', data=lat)
    
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('tclw', data=tclw)
    p.create_dataset('tciw', data=tciw)
    
    p.create_dataset('alt', data=alt)
    p.create_dataset('alt_c', data=alt_c)
    p.create_dataset('alt_t', data=alt_t)

    p.create_dataset('air_density_g', data=air_density_g)
    p.create_dataset('air_density_so', data=air_density_so)

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

    p.create_dataset('cf_t', data=cf_t)
    p.create_dataset('cf_t_so', data=cf_t_so)
    p.create_dataset('lw_t', data=lw_t)
    p.create_dataset('lw_t_so', data=lw_t_so)
    p.create_dataset('iw_t', data=iw_t)
    p.create_dataset('iw_t_so', data=iw_t_so)

    p.close()






















