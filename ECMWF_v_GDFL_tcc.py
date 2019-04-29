# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 20:13:24 2019

@author: Tristan O'Hanlon

Run:
ECMWF/ECMWF_lat_phase.py
GFDL/GFDL_lat_phase.py

To ensure that the following variables are in memory:
ecmwf_tcc_lat    
gfdl_tcc_lat

[:,0] is latitude
[:,1] is total cloud cover
"""
import matplotlib.pyplot as plt

plt.figure()

fig, ax = plt.subplots()
ax.plot(ecmwf_tcc_lat[:,0],ecmwf_tcc_lat[:,1], '-r', label='ECMWF')
ax.plot(gfdl_tcc_lat[:,0],gfdl_tcc_lat[:,1], '-g', label='GFDL')
#ax.axis('equal')
leg = ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

plt.xlabel('Latitude')
plt.ylabel('Cloud Cover')
plt.title('Total Cloud Fraction vs Latitude - 1979 - Present')
plt.grid(True)
plt.show()





















