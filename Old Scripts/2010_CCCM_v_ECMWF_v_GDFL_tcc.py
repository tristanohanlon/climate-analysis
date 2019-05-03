# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 20:13:24 2019

@author: Tristan O'Hanlon

Run:
ECMWF/ECMWF_lat_phase.py - ensure that only 2010 values are selected
GFDL/GFDL_lat_phase.py - ensure that only 2010 values are selected
CERES/2010_CERES_lat_tcc.py
CCCM/2010_CCCM_lat_tcc.py

To ensure that the following variables are in memory:
ecmwf_tcc_lat    
gfdl_tcc_lat
ceres_tcc_lat
cccm_tcc_lat

[:,0] is latitude
[:,1] is total cloud cover
"""

import matplotlib.pyplot as plt

plt.figure()

fig, ax = plt.subplots()
ax.plot(cccm_tcc_lat[:,0],cccm_tcc_lat[:,1], '-b', label='CCCM')
ax.plot(ecmwf_tcc_lat[:,0],ecmwf_tcc_lat[:,1], '-r', label='ECMWF Reanalysis')
ax.plot(gfdl_tcc_lat[:,0],gfdl_tcc_lat[:,1], '-g', label='GFDL.AM4')
#ax.axis('equal')
leg = ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

plt.xlabel('Latitude')
plt.ylabel('Cloud Cover Fraction')
plt.title('Total Cloud Fraction vs Latitude - Southern Ocean - 2010')
plt.grid(True)
plt.show()





















