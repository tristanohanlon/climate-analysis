# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 20:13:24 2019

@author: Tristan
"""

import datetime as dt  # Python standard library datetime  module
import numpy as np
import matplotlib.pyplot as plt

plt.figure()

fig, ax = plt.subplots()
ax.plot(tclw[:,0],tclw[:,1], '-b', label='CERES')
#ax.plot(cccm_tcc_lat[:,0],cccm_tcc_lat[:,1], '-y', label='C3M A-Train')
ax.plot(ecmwf_tclw_lat[:,0],ecmwf_tclw_lat[:,1], '-r', label='ECMWF')
ax.plot(gfdl_tclw_lat[:,0],gfdl_tclw_lat[:,1], '-g', label='GFDL')
#ax.axis('equal')
leg = ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

plt.xlabel('Latitude')
plt.ylabel('Cloud Liquid Water Content Fraction')
plt.title('Cloud Liquid Water Content Fraction vs Latitude - 2010')
plt.grid(True)
plt.show()





















