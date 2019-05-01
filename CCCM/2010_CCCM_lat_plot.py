# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon
"""
import matplotlib.pyplot as plt
import math
import h5py
f = h5py.File('2010_CCCM_global_cloud_phase_lat_data.h5', 'r')

tcc = f['tcc'][:]
tclw = f['tclw'][:] # in g/m^3 
tciw = f['tciw'][:] # in g/m^3
tclw_frac = f['tclw_frac'][:]
tciw_frac = f['tciw_frac'][:]

# Plot cloud fraction, liquid and ice water content vs latitude.
plt.figure()
fig, ax = plt.subplots()

ax.plot(tcc[:,0], tcc[:,1], '-r', label='Total Cloud Cover')
#ax.plot(tclw[:,0], tclw[:,1], '--g', label='Liquid Water Content ($gm^{-3}$)')
#ax.plot(tciw[:,0], tciw[:,1], '--b', label='Ice Water Content ($gm^{-3}$)')
ax.plot(tclw_frac[:,0], tclw_frac[:,1], '--g', label='Liquid Water Fraction')
ax.plot(tciw_frac[:,0], tciw_frac[:,1], '--b', label='Ice Water Content Fraction')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);
           
ax.set_ylabel('Cloud Fraction', color='r')
ax.set_xlabel('Latitude')

plt.title('Global Cloud Cover and Phase vs Latitude - CloudSAT 2010')

plt.grid(True)
plt.show()