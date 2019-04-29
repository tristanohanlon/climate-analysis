# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon
"""
import matplotlib.pyplot as plt
import h5py
p = h5py.File('2010_CCCM_global_cloud_phase_profile_data.h5', 'r')
f = h5py.File('2010_CCCM_global_temp_pressure_profile_data.h5', 'r')

cf = p['Cloud Fraction Profile'][:]
lw = p['Liquid Water Content Profile'][:]
iw = p['Ice Water Content Height Profile'][:]
alt = p['Liquid and Ice Content Height Profile'][:]
alt_c = p['Cloud Height Profile'][:]
alt_t = f['Temperature and Pressure Height Profile'][:]
temp = f['Temperature Profile'][:]
pressure = f['Pressure Profile'][:]


# Plot cloud fraction, liquid and ice water content vs altitude.
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()

ax1.plot(cf[20:113],alt_c[20:113], '-r', label='Cloud Fraction')
ax2.plot(lw[40:135], alt[40:135], '--g', label='Liquid Water')
ax2.plot(iw[40:135], alt[40:135], '--b', label='Ice Water')
#ax.axis('equal')
ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);
           
ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
         ncol=4, fancybox=True, shadow=True);
           
ax1.set_xlabel('Cloud Fraction (%)', color='r')
ax1.set_ylabel('Altitude (km)')
ax2.set_xlabel('Liquid and Ice Water Content ($gm^{-3}$)', color='b')

plt.title('Global Cloud Fraction, Ice and Water Content vs Altitude - CALIPSO-CloudSAT 2010')

plt.grid(True)
plt.show()

"""
# Plot temperature and pressure vs altitude.
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()

ax1.plot(temp[1:130],alt_t[1:130], '-r', label='Temperature')
ax2.plot(pressure[1:130], alt_t[1:130], '-b', label='Pressure')
#ax.axis('equal')
ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);
           
ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
         ncol=4, fancybox=True, shadow=True);
           
ax1.set_xlabel('Temperature (K)', color='r')
ax1.set_ylabel('Altitude (km)')
ax2.set_xlabel('Pressure (hPa)', color='b')

plt.title('Global Temperature and Pressure vs Altitude - CALIPSO-CloudSAT 2010')

plt.grid(True)
plt.show()
"""