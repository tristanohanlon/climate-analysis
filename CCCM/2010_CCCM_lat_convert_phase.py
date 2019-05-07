# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon
"""
import h5py
import os
import numpy as np


os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
f = h5py.File('2010_CCCM_so_profile.h5', 'r')

cf = f['Southern Ocean Cloud Fraction Profile'][:]
lw = f['Southern Ocean Liquid Water Content Profile'][:]
iw = f['Southern Ocean Ice Water Content Profile'][:]
temp = f['Southern Ocean Temperature Profile'][:]
pressure = f['Southern Ocean Pressure Profile'][:]

iw = iw[1:110]

air_density = [] #create empty list
#calculate air density at each pressure layer
air_density = (pressure[:,1] * 100) / (286.9 * temp[:,1]) * 1000
air_density = air_density[26:135]

lwc = lw[:,1] / air_density
iwc = iw[:,1] / air_density

temp_cf = temp[35:135]
templiwc = temp[26:135]

cf_temp = np.vstack((temp_cf[:,1],cf[:,1])).T
lwc_temp = np.vstack((templiwc[:,1],lwc)).T
iwc_temp = np.vstack((templiwc[:,1],iwc)).T

pres_cf = pressure[35:135]
presliwc = pressure[26:135]

cf_pres = np.vstack((pres_cf[:,1],cf[:,1])).T
lwc_pres = np.vstack((presliwc[:,1],lwc)).T
iwc_pres = np.vstack((presliwc[:,1],lwc)).T


os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
f = h5py.File('2010_CCCM_global_latitude.h5', 'r')

tcc = f['Cloud Fraction'][:]
tclw = f['Liquid Water Content'][:] # in g/m^3 
tciw = f['Ice Water Content'][:] # in g/m^3

# Plot cloud fraction, liquid and ice water content vs latitude.
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(tcc[:,0], tcc[:,1], '-r', label='Total Cloud Cover')
ax2.plot(tclw[:,0], tclw[:,1], '-b', label='Liquid Water Content ($gm^{-3}$)')
ax2.plot(tciw[:,0], tciw[:,1], '--b', label='Ice Water Content ($gm^{-3}$)')

ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);
ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
          ncol=4, fancybox=True, shadow=True);  
ax1.set_ylabel('Cloud Fraction', color='r')
ax1.set_xlabel('Latitude')
ax2.set_ylabel('Liquid and Ice Water Content ($gm^{-3}$)', color='r')

plt.title('Global Cloud Cover and Phase Content vs Latitude - CCCM 2010')

plt.grid(True)
plt.show()