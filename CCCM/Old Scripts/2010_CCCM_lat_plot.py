# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon
"""
import matplotlib.pyplot as plt
import math
import h5py
import os
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
f = h5py.File('2010_CCCM_global_latitude.h5', 'r')

tcc = f['Cloud Fraction'][:]
tclw = f['Liquid Water Content'][:] # in g/m^3 
tciw = f['Ice Water Content'][:] # in g/m^3

# Plot cloud fraction, liquid and ice water content vs latitude.
plt.figure()
fig, ax = plt.subplots()

#ax2 = ax1.twinx()
ax.plot(cccm_tclw_lat[:,0], cccm_tclw_lat[:,1], '-r', label='Liquid Water Content')
ax.plot(cccm_tciw_lat[:,0], cccm_tciw_lat[:,1], '--r', label='Ice Water Content')
ax.plot(cccm_tclw_av_lat[:,0], cccm_tclw_av_lat[:,1], '-b', label='Mean Liquid Water Content')
ax.plot(cccm_tciw_av_lat[:,0], cccm_tciw_av_lat[:,1], '--b', label='Mean Ice Water Content')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);
#ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
#          ncol=4, fancybox=True, shadow=True);  
ax.set_ylabel('Liquid and Ice Water Content', color='r')
ax.set_xlabel('Latitude')
#ax2.set_ylabel('Liquid and Ice Water Content ($gm^{-3}$)', color='r')

plt.title('Global Cloud Cover and Phase Content vs Latitude - CCCM 2010')

plt.grid(True)
plt.show()