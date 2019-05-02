# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon
"""
import matplotlib.pyplot as plt
import math
import h5py
import os
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM')
f = h5py.File('2010_CCCM_SO_profile_data_test1.h5', 'r')

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