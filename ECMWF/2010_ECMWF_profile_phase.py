# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

This will create a dataset of time averaged global cloud cover, liquid water cloud fraction and ice water cloud fraction with altitude.
Data is stored in the 2D arrays: 

ecmwf_tcc_alt
ecmwf_tclw_alt
ecmwf_tciw_alt

[:,0] = alt
[:,1] = cloud fraction
"""

import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.pyplot as plt


dataset = Dataset('Data/pressure_levels/2010_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
#lon = dataset.variables['longitude'][:] #Extract longitude data
plevel = dataset.variables['level'][:] #Extract pressure (millibars) level data (37 levels)
tcc = dataset.variables['cc'][:] #Extract fraction of cloud cover, keyed to time, lat and lon (12, 37, 721, 1440)
tciw = dataset.variables['ciwc'][:] #Extract specific cloud ice water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
tclw = dataset.variables['clwc'][:] #Extract specific cloud liquid water content (kg/kg), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)
temp = dataset.variables['t'][:] #Extract air temperature profile content (K), keyed to time, pressure level, lat and lon (12, 37, 721, 1440)

tcc = np.array(tcc)
tcc = np.mean(tcc, axis=0) #Average fraction of cloud cover over time
tcc = np.mean(tcc, axis=-1) #Average fraction of cloud cover over longitude
tcc = np.mean(tcc, axis=-1) #Average fraction of cloud cover over latitude

tciw = np.array(tciw)
tciw = np.mean(tciw, axis=0) #Average specific cloud ice water content (kg/kg) over time
tciw = np.mean(tciw, axis=-1) #Average specific cloud ice water content (kg/kg) over longitude
tciw = np.mean(tciw, axis=-1) * 10000 #Average specific cloud ice water content (kg/kg) over latitude, scale up by 10000

tclw = np.array(tclw)
tclw = np.mean(tclw, axis=0) #Average specific cloud liquid water content (kg/kg) over time
tclw = np.mean(tclw, axis=-1) #Average specific cloud liquid water content (kg/kg) over longitude
tclw = np.mean(tclw, axis=-1) * 10000 #Average specific cloud liquid water content (kg/kg) over latitude, scale up by 10000

temp = np.array(temp)
temp = np.mean(temp, axis=0) #Average air temperature over time
temp = np.mean(temp, axis=-1) #Average air temperature over longitude
temp = np.mean(temp, axis=-1)  #Average air temperature over latitude

# Join the two lists as if they were two columns side by side, into a list of two elements each
ecmwf_tcc_plevel = np.vstack((plevel, tcc)).T 
ecmwf_tclw_plevel = np.vstack((plevel, tclw)).T
ecmwf_tciw_plevel = np.vstack((plevel, tciw)).T
ecmwf_temp_plevel = np.vstack((plevel, temp)).T

#----------------------------#

plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()
ax1.plot(ecmwf_tcc_plevel[:,1],ecmwf_tcc_plevel[:,0], '-r', label='Fraction Cloud Cover')
ax2.plot(ecmwf_tclw_plevel[:,1],ecmwf_tclw_plevel[:,0], '-b', label='Specific Cloud Liquid Water Content')
ax2.plot(ecmwf_tciw_plevel[:,1],ecmwf_tciw_plevel[:,0], '--b', label='Specific Cloud Ice Water Content')

#ax.axis('equal')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Specific Cloud Liquid and Ice Water Content (kg/kg) x $10^{-4}$')
ax1.set_ylabel('Pressure Level (hPa)')
plt.title('Cloud Fraction and Phase vs Pressure Level ECMWF 1979 - Present')
plt.gca().invert_yaxis()

plt.grid(True)
plt.show()