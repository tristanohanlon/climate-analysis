# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

Takes the reduced and averaged datasets and plots the profile and latitude data.
"""
import matplotlib.pyplot as plt
import h5py
p = h5py.File('2010_CERES_global_cloud_phase_profile_data.h5', 'r')
f = h5py.File('2010_CERES_global_cloud_phase_lat_data.h5', 'r')

cfp = p['Cloud Fraction Profile'][:]
lwp = p['Liquid Water Content Profile'][:]
iwp = p['Ice Water Content Height Profile'][:]

tcc = f['Total Cloud Cover'][:]
tciw = f['Total Ice Water Cloud Cover'][:]
tclw = f['Total Liquid Water Cloud Cover'][:]

# Plot cloud fraction, liquid and ice water content vs altitude.
plt.figure()
fig, ax = plt.subplots()


ax.plot(cfp[:,1],cfp[:,0], '-r', label='Cloud Fraction')
ax.plot(lwp[:,1],lwp[:,0], '--g', label='Liquid Water')
ax.plot(iwp[:,1],iwp[:,0], '--b', label='Ice Water')
#ax.axis('equal')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

           
ax.set_xlabel('Cloud Fraction', color='r')
ax.set_ylabel('Altitude (km)')

plt.title('Global Cloud Fraction and Phase vs Altitude - CERES 2010')

plt.grid(True)
plt.show()


# Plot cloud cover and phase vs latitude.
plt.figure()
fig, ax = plt.subplots()


ax.plot(tcc[:,0],tcc[:,1], '-r', label='Cloud Fraction')
ax.plot(tclw[:,0],tclw[:,1], '--g', label='Liquid Water')
ax.plot(tciw[:,0],tciw[:,1], '--b', label='Ice Water')
#ax.axis('equal')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

           
ax.set_ylabel('Cloud Fraction', color='r')
ax.set_xlabel('Latitude')

plt.title('Global Cloud Fraction and Phase vs Latitude - CERES 2010')

plt.grid(True)
plt.show()
