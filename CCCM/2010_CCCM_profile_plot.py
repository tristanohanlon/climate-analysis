# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon
"""
import matplotlib.pyplot as plt
import h5py
import os
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM')
#f = h5py.File('2010_CCCM_SO_lat_alt_profile_data.h5', 'r')
f = h5py.File('2010_CCCM_global_cloud_phase_temp_pressure_profile.h5', 'r')



#lat = f['Latitude'][:]
#alt = f['Phase Altitude'][:]
#alt_t = f['Pressure and Temperature Altitude'][:]

#cf = f['Southern Ocean Cloud Fraction Profile'][:]
#lw = f['Southern Ocean Liquid Water Content Profile'][:]
#iw = f['Southern Ocean Ice Water Content Profile'][:]
#temp = f['Southern Ocean Temperature Profile'][:]
#pressure = f['Southern Ocean Pressure Profile'][:]

cf = f['Cloud Fraction Profile'][:]
lw = f['Liquid Water Content Profile'][:]
iw = f['Ice Water Content Profile'][:]
temp = f['Temperature Profile'][:]
pressure = f['Pressure Profile'][:]


#Select only data between 0 and 20km
#cf = cf[cf[:,0]<=20]
#cf = cf[cf[:,0]>=0]
#lw = lw[lw[:,0]<=20]
#lw = lw[lw[:,0]>=0]
#iw = iw[iw[:,0]<=20]
#iw = iw[iw[:,0]>=0]
"""
# Plot cloud fraction, liquid and ice water content vs altitude.
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()

ax1.plot(cf[:,1],cf[:,0], '-r', label='Cloud Fraction')
ax2.plot(lw[:,1], lw[:,0], '--g', label='Liquid Water Content')
ax2.plot(iw[:,1], iw[:,0], '--b', label='Ice Water Content')
#ax.axis('equal')
ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);
           
ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
         ncol=4, fancybox=True, shadow=True);
           
ax1.set_xlabel('Cloud Fraction (%)', color='r')
ax1.set_ylabel('Altitude (km)')
ax2.set_xlabel('Liquid and Ice Water Content ($gm^{-3}$)', color='b')

plt.title('Southern Ocean Cloud Fraction, Ice and Water Content vs Altitude - CALIPSO-CloudSAT 2010')

plt.grid(True)
plt.show()


# Plot temperature and pressure vs altitude.
plt.figure()
fig, ax1 = plt.subplots()

#ax2 = ax1.twiny()

ax1.plot(temp[:,1],temp[:,0], '-r', label='Temperature')
#ax2.plot(pressure[:,1], pressure[:,0], '-b', label='Pressure')
#ax.axis('equal')
ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);
           
#ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
#         ncol=4, fancybox=True, shadow=True);
           
ax1.set_xlabel('Temperature (K)', color='r')
ax1.set_ylabel('Altitude (km)')
#ax2.set_xlabel('Pressure (hPa)', color='b')

plt.title('Southern Ocean Temperature vs Altitude - CALIPSO-CloudSAT 2010')

plt.grid(True)
plt.show()
"""