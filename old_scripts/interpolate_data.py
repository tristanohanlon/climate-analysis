# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 17:30:26 2019

@author: Tristan
"""
import time
import sys
import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import h5py
import math
from scipy import integrate
from scipy import interpolate



#-----------Latitude-----------#
lat = cccm_tcc_lat_g_enhanced[:,0]


f = interpolate.interp1d(ceres_tcc_lat_g[:,0], ceres_tcc_lat_g[:,1], fill_value="extrapolate", kind = 'cubic')
ceres_tcc_lat = f(lat)
ceres_tcc_lat_g = np.vstack((lat, ceres_tcc_lat)).T

f = interpolate.interp1d(ecmwf_tcc_lat_g[:,0], ecmwf_tcc_lat_g[:,1], fill_value="extrapolate", kind = 'cubic')
ecmwf_tcc_lat = f(lat)
ecmwf_tcc_lat_g = np.vstack((lat, ecmwf_tcc_lat)).T

f = interpolate.interp1d(gfdl4_tcc_lat_g[:,0], gfdl4_tcc_lat_g[:,1], fill_value="extrapolate", kind = 'cubic')
gfdl4_tcc_lat = f(lat)
gfdl4_tcc_lat_g = np.vstack((lat, gfdl4_tcc_lat)).T

f = interpolate.interp1d(gfdl_hiram_tcc_lat_g[:,0], gfdl_hiram_tcc_lat_g[:,1], fill_value="extrapolate", kind = 'cubic')
gfdl_hiram_tcc_lat = f(lat)
gfdl_hiram_tcc_lat_g = np.vstack((lat, gfdl_hiram_tcc_lat)).T

f = interpolate.interp1d(cam5_tcc_lat_g[:,0], cam5_tcc_lat_g[:,1], fill_value="extrapolate", kind = 'cubic')
cam5_tcc_lat = f(lat)
cam5_tcc_lat_g = np.vstack((lat, cam5_tcc_lat)).T

f = interpolate.interp1d(cam6_tcc_lat_g[:,0], cam6_tcc_lat_g[:,1], fill_value="extrapolate", kind = 'cubic')
cam6_tcc_lat = f(lat)
cam6_tcc_lat_g = np.vstack((lat, cam6_tcc_lat)).T

f = interpolate.interp1d(mri_cgcm_tcc_lat_g[:,0], mri_cgcm_tcc_lat_g[:,1], fill_value="extrapolate", kind = 'cubic')
mri_cgcm_tcc_lat = f(lat)
mri_cgcm_tcc_lat_g = np.vstack((lat, mri_cgcm_tcc_lat)).T

f = interpolate.interp1d(mri_tcc_lat_g[:,0], mri_tcc_lat_g[:,1], fill_value="extrapolate", kind = 'cubic')
mri_tcc_lat = f(lat)
mri_tcc_lat_g = np.vstack((lat, mri_tcc_lat)).T

f = interpolate.interp1d(miroc5_tcc_lat_g[:,0], miroc5_tcc_lat_g[:,1], fill_value="extrapolate", kind = 'cubic')
miroc5_tcc_lat = f(lat)
miroc5_tcc_lat_g = np.vstack((lat, miroc5_tcc_lat)).T

f = interpolate.interp1d(miroc6_tcc_lat_g[:,0], miroc6_tcc_lat_g[:,1], fill_value="extrapolate", kind = 'cubic')
miroc6_tcc_lat = f(lat)
miroc6_tcc_lat_g = np.vstack((lat, miroc6_tcc_lat)).T

f = interpolate.interp1d(ipsl5_tcc_lat_g[:,0], ipsl5_tcc_lat_g[:,1], fill_value="extrapolate", kind = 'cubic')
ipsl5_tcc_lat = f(lat)
ipsl5_tcc_lat_g = np.vstack((lat, ipsl5_tcc_lat)).T

f = interpolate.interp1d(ipsl6_tcc_lat_g[:,0], ipsl6_tcc_lat_g[:,1], fill_value="extrapolate", kind = 'cubic')
ipsl6_tcc_lat = f(lat)
ipsl6_tcc_lat_g = np.vstack((lat, ipsl6_tcc_lat)).T

f = interpolate.interp1d(giss5_tcc_lat_g[:,0], giss5_tcc_lat_g[:,1], fill_value="extrapolate", kind = 'cubic')
giss5_tcc_lat = f(lat)
giss5_tcc_lat_g = np.vstack((lat, giss5_tcc_lat)).T

f = interpolate.interp1d(giss6_tcc_lat_g[:,0], giss6_tcc_lat_g[:,1], fill_value="extrapolate", kind = 'cubic')
giss6_tcc_lat = f(lat)
giss6_tcc_lat_g = np.vstack((lat, giss6_tcc_lat)).T


cmip5_lat_max = np.maximum.reduce([ecmwf_tcc_lat, gfdl_hiram_tcc_lat, cam5_tcc_lat, mri_cgcm_tcc_lat, miroc5_tcc_lat, ipsl5_tcc_lat, giss5_tcc_lat])
cmip5_lat_min = np.minimum.reduce([ecmwf_tcc_lat, gfdl_hiram_tcc_lat, cam5_tcc_lat, mri_cgcm_tcc_lat, miroc5_tcc_lat, ipsl5_tcc_lat, giss5_tcc_lat])

cmip6_lat_max = np.maximum.reduce([ecmwf_tcc_lat, gfdl4_tcc_lat, cam6_tcc_lat, mri_tcc_lat, miroc6_tcc_lat, ipsl6_tcc_lat, giss6_tcc_lat])
cmip6_lat_min = np.minimum.reduce([ecmwf_tcc_lat, gfdl4_tcc_lat, cam6_tcc_lat, mri_tcc_lat, miroc6_tcc_lat, ipsl6_tcc_lat, giss6_tcc_lat])



"""
fig, ax1 = plt.subplots()

ax1.plot(ceres_tcc_lat_g[:,0],ceres_tcc_lat_g[:,1], '-b', label='CERES - Satellite', linewidth=1)
ax1.plot(gfdl_hiram_tcc_lat_g[:,0],gfdl_hiram_tcc_lat_g[:,1], '--g', label='CMIP5-GFDL-HIRAM-AMIP', linewidth=1)
ax1.plot(cam5_tcc_lat_g[:,0],cam5_tcc_lat_g[:,1], '--m', label='CMIP5-CESM1-CAM5-AMIP', linewidth=1)


ax1.plot(gfdl4_tcc_lat_g[:,0],gfdl4_tcc_lat_g[:,1], '-g', label='CMIP6-GFDL-AM4-AMIP', linewidth=1)
ax1.plot(cam6_tcc_lat_g[:,0],cam6_tcc_lat_g[:,1], '-m', label='CMIP6-CESM2-CAM6-AMIP', linewidth=1)


ax1.fill_between(lat, cmip5_lat_min, cmip5_lat_max, facecolor='black', alpha=0.2,
                label='CMIP5 Model Range')

ax1.fill_between(lat, cmip6_lat_min, cmip6_lat_max, facecolor='yellow', alpha=0.4,
                label='CMIP6 Model Range')


ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));


ax1.set_ylabel('Cloud Fraction')
ax1.set_xlabel('Latitude')

#ax1.set_title ('Global Cloud Fraction vs Latitude')

ax1.grid(True)
plt.savefig("tcc_lat_g.pdf", format="pdf", bbox_inches='tight')
plt.show()
"""

#-----------TCC_g altitude-----------#

alt = cccm_tcc_alt_g[7:100,0]


f = interpolate.interp1d(calipso_tcc_alt_g[:,0], calipso_tcc_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
calipso_tcc_alt = f(alt)
calipso_tcc_alt_g = np.vstack((alt, calipso_tcc_alt)).T

f = interpolate.interp1d(ecmwf_tcc_alt_g[:,0], ecmwf_tcc_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
ecmwf_tcc_alt = f(alt)
ecmwf_tcc_alt_g = np.vstack((alt, ecmwf_tcc_alt)).T

f = interpolate.interp1d(gfdl4_tcc_alt_g[:,0], gfdl4_tcc_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
gfdl4_tcc_alt = f(alt)
gfdl4_tcc_alt_g = np.vstack((alt, gfdl4_tcc_alt)).T

f = interpolate.interp1d(gfdl_hiram_tcc_alt_g[:,0], gfdl_hiram_tcc_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
gfdl_hiram_tcc_alt = f(alt)
gfdl_hiram_tcc_alt_g = np.vstack((alt, gfdl_hiram_tcc_alt)).T

f = interpolate.interp1d(cam5_tcc_alt_g[:,0], cam5_tcc_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
cam5_tcc_alt = f(alt)
cam5_tcc_alt_g = np.vstack((alt, cam5_tcc_alt)).T

f = interpolate.interp1d(cam6_tcc_alt_g[:,0], cam6_tcc_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
cam6_tcc_alt = f(alt)
cam6_tcc_alt_g = np.vstack((alt, cam6_tcc_alt)).T

f = interpolate.interp1d(mri_cgcm_tcc_alt_g[:,0], mri_cgcm_tcc_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
mri_cgcm_tcc_alt = f(alt)
mri_cgcm_tcc_alt_g = np.vstack((alt, mri_cgcm_tcc_alt)).T

f = interpolate.interp1d(mri_tcc_alt_g[:,0], mri_tcc_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
mri_tcc_alt = f(alt)
mri_tcc_alt_g = np.vstack((alt, mri_tcc_alt)).T

f = interpolate.interp1d(miroc5_tcc_alt_g[:,0], miroc5_tcc_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
miroc5_tcc_alt = f(alt)
miroc5_tcc_alt_g = np.vstack((alt, miroc5_tcc_alt)).T

f = interpolate.interp1d(miroc6_tcc_alt_g[:,0], miroc6_tcc_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
miroc6_tcc_alt = f(alt)
miroc6_tcc_alt_g = np.vstack((alt, miroc6_tcc_alt)).T

f = interpolate.interp1d(ipsl5_tcc_alt_g[:,0], ipsl5_tcc_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
ipsl5_tcc_alt = f(alt)
ipsl5_tcc_alt_g = np.vstack((alt, ipsl5_tcc_alt)).T

f = interpolate.interp1d(ipsl6_tcc_alt_g[:,0], ipsl6_tcc_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
ipsl6_tcc_alt = f(alt)
ipsl6_tcc_alt_g = np.vstack((alt, ipsl6_tcc_alt)).T

f = interpolate.interp1d(giss5_tcc_alt_g[:,0], giss5_tcc_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
giss5_tcc_alt = f(alt)
giss5_tcc_alt_g = np.vstack((alt, giss5_tcc_alt)).T

f = interpolate.interp1d(giss6_tcc_alt_g[:,0], giss6_tcc_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
giss6_tcc_alt = f(alt)
giss6_tcc_alt_g = np.vstack((alt, giss6_tcc_alt)).T


cmip5_alt_max = np.maximum.reduce([ecmwf_tcc_alt, gfdl_hiram_tcc_alt, cam5_tcc_alt, mri_cgcm_tcc_alt, miroc5_tcc_alt, ipsl5_tcc_alt, giss5_tcc_alt])
cmip5_alt_min = np.minimum.reduce([ecmwf_tcc_alt, gfdl_hiram_tcc_alt, cam5_tcc_alt, mri_cgcm_tcc_alt, miroc5_tcc_alt, ipsl5_tcc_alt, giss5_tcc_alt])

cmip6_alt_max = np.maximum.reduce([ecmwf_tcc_alt, gfdl4_tcc_alt, cam6_tcc_alt, mri_tcc_alt, miroc6_tcc_alt, ipsl6_tcc_alt, giss6_tcc_alt])
cmip6_alt_min = np.minimum.reduce([ecmwf_tcc_alt, gfdl4_tcc_alt, cam6_tcc_alt, mri_tcc_alt, miroc6_tcc_alt, ipsl6_tcc_alt, giss6_tcc_alt])


"""
fig, ax1 = plt.subplots()

ax1.plot(cccm_tcc_alt_g[7:100,1],cccm_tcc_alt_g[7:100,0], '-b', label='CCCM - Satellite', linewidth=1)

ax1.plot(gfdl_hiram_tcc_alt_g[:,1],gfdl_hiram_tcc_alt_g[:,0], '--g', label='CMIP5-GFDL-HIRAM-AMIP', linewidth=1)
ax1.plot(cam5_tcc_alt_g[:,1],cam5_tcc_alt_g[:,0], '--m', label='CMIP5-CESM1-CAM5-AMIP', linewidth=1)

ax1.plot(gfdl4_tcc_alt_g[:,1],gfdl4_tcc_alt_g[:,0], '-g', label='CMIP6-GFDL-AM4-AMIP', linewidth=1)
ax1.plot(cam6_tcc_alt_g[:,1],cam6_tcc_alt_g[:,0], '-m', label='CMIP6-CESM2-CAM6-AMIP', linewidth=1)

ax1.fill_betweenx(alt, cmip5_alt_min, cmip5_alt_max, facecolor='black', alpha=0.2,
                label='CMIP5 Model Range')

ax1.fill_betweenx(alt, cmip6_alt_min, cmip6_alt_max, facecolor='yellow', alpha=0.4,
                label='CMIP6 Model Range')



ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Fraction')

#ax1.set_title ('Global Cloud Fraction vs Latitude')
ax1.set_xlim(0, 0.6)
ax1.grid(True)
plt.savefig("tcc_alt_g.pdf", format="pdf", bbox_inches='tight')
plt.show()
"""

#-----------TCC_so altitude-----------#

alt_so = cccm_tcc_alt_so[7:100,0]


f = interpolate.interp1d(calipso_tcc_alt_so[:,0], calipso_tcc_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
calipso_tcc_alt_s = f(alt_so)
calipso_tcc_alt_so = np.vstack((alt_so, calipso_tcc_alt_s)).T

f = interpolate.interp1d(ecmwf_tcc_alt_so[:,0], ecmwf_tcc_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
ecmwf_tcc_alt_s = f(alt_so)
ecmwf_tcc_alt_so = np.vstack((alt_so, ecmwf_tcc_alt_s)).T

f = interpolate.interp1d(gfdl4_tcc_alt_so[:,0], gfdl4_tcc_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
gfdl4_tcc_alt_s = f(alt_so)
gfdl4_tcc_alt_so = np.vstack((alt_so, gfdl4_tcc_alt_s)).T

f = interpolate.interp1d(gfdl_hiram_tcc_alt_so[:,0], gfdl_hiram_tcc_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
gfdl_hiram_tcc_alt_s = f(alt_so)
gfdl_hiram_tcc_alt_so = np.vstack((alt_so, gfdl_hiram_tcc_alt_s)).T

f = interpolate.interp1d(cam5_tcc_alt_so[:,0], cam5_tcc_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
cam5_tcc_alt_s = f(alt_so)
cam5_tcc_alt_so = np.vstack((alt_so, cam5_tcc_alt_s)).T

f = interpolate.interp1d(cam6_tcc_alt_so[:,0], cam6_tcc_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
cam6_tcc_alt_s = f(alt_so)
cam6_tcc_alt_so = np.vstack((alt_so, cam6_tcc_alt_s)).T

f = interpolate.interp1d(mri_cgcm_tcc_alt_so[:,0], mri_cgcm_tcc_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
mri_cgcm_tcc_alt_s = f(alt_so)
mri_cgcm_tcc_alt_so = np.vstack((alt_so, mri_cgcm_tcc_alt_s)).T

f = interpolate.interp1d(mri_tcc_alt_so[:,0], mri_tcc_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
mri_tcc_alt_s = f(alt_so)
mri_tcc_alt_so = np.vstack((alt_so, mri_tcc_alt_s)).T

f = interpolate.interp1d(miroc5_tcc_alt_so[:,0], miroc5_tcc_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
miroc5_tcc_alt_s = f(alt_so)
miroc5_tcc_alt_so = np.vstack((alt_so, miroc5_tcc_alt_s)).T

f = interpolate.interp1d(miroc6_tcc_alt_so[:,0], miroc6_tcc_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
miroc6_tcc_alt_s = f(alt_so)
miroc6_tcc_alt_so = np.vstack((alt_so, miroc6_tcc_alt_s)).T

f = interpolate.interp1d(ipsl5_tcc_alt_so[:,0], ipsl5_tcc_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
ipsl5_tcc_alt_s = f(alt_so)
ipsl5_tcc_alt_so = np.vstack((alt_so, ipsl5_tcc_alt_s)).T

f = interpolate.interp1d(ipsl6_tcc_alt_so[:,0], ipsl6_tcc_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
ipsl6_tcc_alt_s = f(alt_so)
ipsl6_tcc_alt_so = np.vstack((alt_so, ipsl6_tcc_alt_s)).T

f = interpolate.interp1d(giss5_tcc_alt_so[:,0], giss5_tcc_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
giss5_tcc_alt_s = f(alt_so)
giss5_tcc_alt_so = np.vstack((alt_so, giss5_tcc_alt_s)).T

f = interpolate.interp1d(giss6_tcc_alt_so[:,0], giss6_tcc_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
giss6_tcc_alt_s = f(alt_so)
giss6_tcc_alt_so = np.vstack((alt_so, giss6_tcc_alt_s)).T


cmip5_alt_max = np.maximum.reduce([ecmwf_tcc_alt_s, gfdl_hiram_tcc_alt_s, cam5_tcc_alt_s, mri_cgcm_tcc_alt_s, miroc5_tcc_alt_s, ipsl5_tcc_alt_s, giss5_tcc_alt_s])
cmip5_alt_min = np.minimum.reduce([ecmwf_tcc_alt_s, gfdl_hiram_tcc_alt_s, cam5_tcc_alt_s, mri_cgcm_tcc_alt_s, miroc5_tcc_alt_s, ipsl5_tcc_alt_s, giss5_tcc_alt_s])

cmip6_alt_max = np.maximum.reduce([ecmwf_tcc_alt_s, gfdl4_tcc_alt_s, cam6_tcc_alt_s, mri_tcc_alt_s, miroc6_tcc_alt_s, ipsl6_tcc_alt_s, giss6_tcc_alt_s])
cmip6_alt_min = np.minimum.reduce([ecmwf_tcc_alt_s, gfdl4_tcc_alt_s, cam6_tcc_alt_s, mri_tcc_alt_s, miroc6_tcc_alt_s, ipsl6_tcc_alt_s, giss6_tcc_alt_s])


"""
fig, ax1 = plt.subplots()

ax1.plot(cccm_tcc_alt_so[7:100,1],cccm_tcc_alt_so[7:100,0], '-b', label='CCCM - Satellite', linewidth=1)

ax1.plot(gfdl_hiram_tcc_alt_so[:,1],gfdl_hiram_tcc_alt_so[:,0], '--g', label='CMIP5-GFDL-HIRAM-AMIP', linewidth=1)
ax1.plot(cam5_tcc_alt_so[:,1],cam5_tcc_alt_so[:,0], '--m', label='CMIP5-CESM1-CAM5-AMIP', linewidth=1)

ax1.plot(gfdl4_tcc_alt_so[:,1],gfdl4_tcc_alt_so[:,0], '-g', label='CMIP6-GFDL-AM4-AMIP', linewidth=1)
ax1.plot(cam6_tcc_alt_so[:,1],cam6_tcc_alt_so[:,0], '-m', label='CMIP6-CESM2-CAM6-AMIP', linewidth=1)

ax1.fill_betweenx(alt_so, cmip5_alt_min, cmip5_alt_max, facecolor='black', alpha=0.2,
                label='CMIP5 Model Range')

ax1.fill_betweenx(alt_so, cmip6_alt_min, cmip6_alt_max, facecolor='yellow', alpha=0.4,
                label='CMIP6 Model Range')



ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Fraction')

#ax1.set_title ('Global Cloud Fraction vs Latitude')
ax1.set_xlim(0, 0.6)
ax1.grid(True)
plt.savefig("tcc_alt_so.pdf", format="pdf", bbox_inches='tight')
plt.show()
"""


#-----------TCLW_G Frac altitude-----------#

alt_lw = cccm_tclw_alt_g[7:52,0]


f = interpolate.interp1d(calipso_tclw_frac_alt_g[:,0], calipso_tclw_frac_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
calipso_tclw_frac_alt = f(alt_lw)
calipso_tclw_frac_alt_g = np.vstack((alt_lw, calipso_tclw_frac_alt)).T

f = interpolate.interp1d(ecmwf_tclw_frac_alt_g[:,0], ecmwf_tclw_frac_alt_g[:,1], fill_value="extrapolate", kind = 'cubic')
ecmwf_tclw_frac_alt = f(alt_lw)
ecmwf_tclw_frac_alt_g = np.vstack((alt_lw, ecmwf_tclw_frac_alt)).T

f = interpolate.interp1d(gfdl4_tclw_frac_alt_g[:,0], gfdl4_tclw_frac_alt_g[:,1], fill_value="extrapolate")
gfdl4_tclw_frac_alt = f(alt_lw)
gfdl4_tclw_frac_alt_g = np.vstack((alt_lw, gfdl4_tclw_frac_alt)).T

f = interpolate.interp1d(gfdl_hiram_tclw_frac_alt_g[:,0], gfdl_hiram_tclw_frac_alt_g[:,1], fill_value="extrapolate")
gfdl_hiram_tclw_frac_alt = f(alt_lw)
gfdl_hiram_tclw_frac_alt_g = np.vstack((alt_lw, gfdl_hiram_tclw_frac_alt)).T

f = interpolate.interp1d(cam5_tclw_frac_alt_g[:,0], cam5_tclw_frac_alt_g[:,1], fill_value="extrapolate")
cam5_tclw_frac_alt = f(alt_lw)
cam5_tclw_frac_alt_g = np.vstack((alt_lw, cam5_tclw_frac_alt)).T

f = interpolate.interp1d(cam6_tclw_frac_alt_g[:,0], cam6_tclw_frac_alt_g[:,1], fill_value="extrapolate")
cam6_tclw_frac_alt = f(alt_lw)
cam6_tclw_frac_alt_g = np.vstack((alt_lw, cam6_tclw_frac_alt)).T

f = interpolate.interp1d(mri_cgcm_tclw_frac_alt_g[:,0], mri_cgcm_tclw_frac_alt_g[:,1], fill_value="extrapolate")
mri_cgcm_tclw_frac_alt = f(alt_lw)
mri_cgcm_tclw_frac_alt_g = np.vstack((alt_lw, mri_cgcm_tclw_frac_alt)).T

f = interpolate.interp1d(mri_tclw_frac_alt_g[:,0], mri_tclw_frac_alt_g[:,1], fill_value="extrapolate")
mri_tclw_frac_alt = f(alt_lw)
mri_tclw_frac_alt_g = np.vstack((alt_lw, mri_tclw_frac_alt)).T

f = interpolate.interp1d(miroc5_tclw_frac_alt_g[:,0], miroc5_tclw_frac_alt_g[:,1], fill_value="extrapolate")
miroc5_tclw_frac_alt = f(alt_lw)
miroc5_tclw_frac_alt_g = np.vstack((alt_lw, miroc5_tclw_frac_alt)).T

f = interpolate.interp1d(miroc6_tclw_frac_alt_g[:,0], miroc6_tclw_frac_alt_g[:,1], fill_value="extrapolate")
miroc6_tclw_frac_alt = f(alt_lw)
miroc6_tclw_frac_alt_g = np.vstack((alt_lw, miroc6_tclw_frac_alt)).T

f = interpolate.interp1d(ipsl5_tclw_frac_alt_g[:,0], ipsl5_tclw_frac_alt_g[:,1], fill_value="extrapolate")
ipsl5_tclw_frac_alt = f(alt_lw)
ipsl5_tclw_frac_alt_g = np.vstack((alt_lw, ipsl5_tclw_frac_alt)).T

f = interpolate.interp1d(ipsl6_tclw_frac_alt_g[:,0], ipsl6_tclw_frac_alt_g[:,1], fill_value="extrapolate")
ipsl6_tclw_frac_alt = f(alt_lw)
ipsl6_tclw_frac_alt_g = np.vstack((alt_lw, ipsl6_tclw_frac_alt)).T

f = interpolate.interp1d(giss5_tclw_frac_alt_g[:,0], giss5_tclw_frac_alt_g[:,1], fill_value="extrapolate")
giss5_tclw_frac_alt = f(alt_lw)
giss5_tclw_frac_alt_g = np.vstack((alt_lw, giss5_tclw_frac_alt)).T

f = interpolate.interp1d(giss6_tclw_frac_alt_g[:,0], giss6_tclw_frac_alt_g[:,1], fill_value="extrapolate")
giss6_tclw_frac_alt = f(alt_lw)
giss6_tclw_frac_alt_g = np.vstack((alt_lw, giss6_tclw_frac_alt)).T


cmip5_alt_max = np.maximum.reduce([ecmwf_tclw_frac_alt, gfdl_hiram_tclw_frac_alt, cam5_tclw_frac_alt, mri_cgcm_tclw_frac_alt, miroc5_tclw_frac_alt, ipsl5_tclw_frac_alt, giss5_tclw_frac_alt])
cmip5_alt_min = np.minimum.reduce([ecmwf_tclw_frac_alt, gfdl_hiram_tclw_frac_alt, cam5_tclw_frac_alt, mri_cgcm_tclw_frac_alt, miroc5_tclw_frac_alt, ipsl5_tclw_frac_alt, giss5_tclw_frac_alt])

cmip6_alt_max = np.maximum.reduce([ecmwf_tclw_frac_alt, gfdl4_tclw_frac_alt, cam6_tclw_frac_alt, mri_tclw_frac_alt, miroc6_tclw_frac_alt, ipsl6_tclw_frac_alt, giss6_tclw_frac_alt])
cmip6_alt_min = np.minimum.reduce([ecmwf_tclw_frac_alt, gfdl4_tclw_frac_alt, cam6_tclw_frac_alt, mri_tclw_frac_alt, miroc6_tclw_frac_alt, ipsl6_tclw_frac_alt, giss6_tclw_frac_alt])


"""
fig, ax1 = plt.subplots()

ax1.plot(calipso_tclw_frac_alt_g[:,1],calipso_tclw_frac_alt_g[:,0], '-b', label='CALIPSO-GOCCP', linewidth=1)
ax1.plot(cccm_tclw_frac_alt_g[7:52,1],cccm_tclw_frac_alt_g[7:52,0], '--b', label='CCCM - Satellite', linewidth=1)

ax1.plot(gfdl_hiram_tclw_frac_alt_g[:,1],gfdl_hiram_tclw_frac_alt_g[:,0], '--g', label='CMIP5-GFDL-HIRAM-AMIP', linewidth=1)
ax1.plot(cam5_tclw_frac_alt_g[:,1],cam5_tclw_frac_alt_g[:,0], '--m', label='CMIP5-CESM1-CAM5-AMIP', linewidth=1)

ax1.plot(gfdl4_tclw_frac_alt_g[:,1],gfdl4_tclw_frac_alt_g[:,0], '-g', label='CMIP6-GFDL-AM4-AMIP', linewidth=1)
ax1.plot(cam6_tclw_frac_alt_g[:,1],cam6_tclw_frac_alt_g[:,0], '-m', label='CMIP6-CESM2-CAM6-AMIP', linewidth=1)

ax1.fill_betweenx(alt_lw, cmip5_alt_min, cmip5_alt_max, facecolor='black', alpha=0.2,
                label='CMIP5 Model Range')

ax1.fill_betweenx(alt_lw, cmip6_alt_min, cmip6_alt_max, facecolor='yellow', alpha=0.4,
                label='CMIP6 Model Range')



ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Liquid Water Fraction')

#ax1.set_title ('Global Cloud Fraction vs Latitude')
ax1.set_xlim(0, 0.6)
ax1.grid(True)
plt.savefig("tclw_frac_alt_g.pdf", format="pdf", bbox_inches='tight')
plt.show()
"""


#-----------TCLW_SO Frac altitude-----------#

alt_lw_so = cccm_tclw_alt_so[7:52,0]


f = interpolate.interp1d(calipso_tclw_frac_alt_so[:,0], calipso_tclw_frac_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
calipso_tclw_frac_alt_s = f(alt_lw_so)
calipso_tclw_frac_alt_so = np.vstack((alt_lw_so, calipso_tclw_frac_alt)).T

f = interpolate.interp1d(ecmwf_tclw_frac_alt_so[:,0], ecmwf_tclw_frac_alt_so[:,1], fill_value="extrapolate", kind = 'cubic')
ecmwf_tclw_frac_alt_s = f(alt_lw_so)
ecmwf_tclw_frac_alt_so = np.vstack((alt_lw_so, ecmwf_tclw_frac_alt_s)).T

f = interpolate.interp1d(gfdl4_tclw_frac_alt_so[:,0], gfdl4_tclw_frac_alt_so[:,1], fill_value="extrapolate")
gfdl4_tclw_frac_alt_s = f(alt_lw_so)
gfdl4_tclw_frac_alt_so = np.vstack((alt_lw_so, gfdl4_tclw_frac_alt_s)).T

f = interpolate.interp1d(gfdl_hiram_tclw_frac_alt_so[:,0], gfdl_hiram_tclw_frac_alt_so[:,1], fill_value="extrapolate")
gfdl_hiram_tclw_frac_alt_s = f(alt_lw_so)
gfdl_hiram_tclw_frac_alt_so = np.vstack((alt_lw_so, gfdl_hiram_tclw_frac_alt_s)).T

f = interpolate.interp1d(cam5_tclw_frac_alt_so[:,0], cam5_tclw_frac_alt_so[:,1], fill_value="extrapolate")
cam5_tclw_frac_alt_s = f(alt_lw_so)
cam5_tclw_frac_alt_so = np.vstack((alt_lw_so, cam5_tclw_frac_alt_s)).T

f = interpolate.interp1d(cam6_tclw_frac_alt_so[:,0], cam6_tclw_frac_alt_so[:,1], fill_value="extrapolate")
cam6_tclw_frac_alt_s = f(alt_lw_so)
cam6_tclw_frac_alt_so = np.vstack((alt_lw_so, cam6_tclw_frac_alt_s)).T

f = interpolate.interp1d(mri_cgcm_tclw_frac_alt_so[:,0], mri_cgcm_tclw_frac_alt_so[:,1], fill_value="extrapolate")
mri_cgcm_tclw_frac_alt_s = f(alt_lw_so)
mri_cgcm_tclw_frac_alt_so = np.vstack((alt_lw_so, mri_cgcm_tclw_frac_alt_s)).T

f = interpolate.interp1d(mri_tclw_frac_alt_so[:,0], mri_tclw_frac_alt_so[:,1], fill_value="extrapolate")
mri_tclw_frac_alt_s = f(alt_lw_so)
mri_tclw_frac_alt_so = np.vstack((alt_lw_so, mri_tclw_frac_alt_s)).T

f = interpolate.interp1d(miroc5_tclw_frac_alt_so[:,0], miroc5_tclw_frac_alt_so[:,1], fill_value="extrapolate")
miroc5_tclw_frac_alt_s = f(alt_lw_so)
miroc5_tclw_frac_alt_so = np.vstack((alt_lw_so, miroc5_tclw_frac_alt_s)).T

f = interpolate.interp1d(miroc6_tclw_frac_alt_so[:,0], miroc6_tclw_frac_alt_so[:,1], fill_value="extrapolate")
miroc6_tclw_frac_alt_s = f(alt_lw_so)
miroc6_tclw_frac_alt_so = np.vstack((alt_lw_so, miroc6_tclw_frac_alt_s)).T

f = interpolate.interp1d(ipsl5_tclw_frac_alt_so[:,0], ipsl5_tclw_frac_alt_so[:,1], fill_value="extrapolate")
ipsl5_tclw_frac_alt_s = f(alt_lw_so)
ipsl5_tclw_frac_alt_so = np.vstack((alt_lw_so, ipsl5_tclw_frac_alt_s)).T

f = interpolate.interp1d(ipsl6_tclw_frac_alt_so[:,0], ipsl6_tclw_frac_alt_so[:,1], fill_value="extrapolate")
ipsl6_tclw_frac_alt_s = f(alt_lw_so)
ipsl6_tclw_frac_alt_so = np.vstack((alt_lw_so, ipsl6_tclw_frac_alt_s)).T

f = interpolate.interp1d(giss5_tclw_frac_alt_so[:,0], giss5_tclw_frac_alt_so[:,1], fill_value="extrapolate")
giss5_tclw_frac_alt_s = f(alt_lw_so)
giss5_tclw_frac_alt_so = np.vstack((alt_lw_so, giss5_tclw_frac_alt_s)).T

f = interpolate.interp1d(giss6_tclw_frac_alt_so[:,0], giss6_tclw_frac_alt_so[:,1], fill_value="extrapolate")
giss6_tclw_frac_alt_s = f(alt_lw_so)
giss6_tclw_frac_alt_so = np.vstack((alt_lw_so, giss6_tclw_frac_alt_s)).T


cmip5_alt_max = np.maximum.reduce([ecmwf_tclw_frac_alt_s, gfdl_hiram_tclw_frac_alt_s, cam5_tclw_frac_alt_s, mri_cgcm_tclw_frac_alt_s, miroc5_tclw_frac_alt_s, ipsl5_tclw_frac_alt_s, giss5_tclw_frac_alt_s])
cmip5_alt_min = np.minimum.reduce([ecmwf_tclw_frac_alt_s, gfdl_hiram_tclw_frac_alt_s, cam5_tclw_frac_alt_s, mri_cgcm_tclw_frac_alt_s, miroc5_tclw_frac_alt_s, ipsl5_tclw_frac_alt_s, giss5_tclw_frac_alt_s])

cmip6_alt_max = np.maximum.reduce([ecmwf_tclw_frac_alt_s, gfdl4_tclw_frac_alt_s, cam6_tclw_frac_alt_s, mri_tclw_frac_alt_s, miroc6_tclw_frac_alt_s, ipsl6_tclw_frac_alt_s, giss6_tclw_frac_alt_s])
cmip6_alt_min = np.minimum.reduce([ecmwf_tclw_frac_alt_s, gfdl4_tclw_frac_alt_s, cam6_tclw_frac_alt_s, mri_tclw_frac_alt_s, miroc6_tclw_frac_alt_s, ipsl6_tclw_frac_alt_s, giss6_tclw_frac_alt_s])


"""
fig, ax1 = plt.subplots()

ax1.plot(calipso_tclw_frac_alt_so[:,1],calipso_tclw_frac_alt_so[:,0], '-b', label='CALIPSO-GOCCP', linewidth=1)
ax1.plot(cccm_tclw_frac_alt_so[7:52,1],cccm_tclw_frac_alt_so[7:52,0], '--b', label='CCCM - Satellite', linewidth=1)

ax1.plot(gfdl_hiram_tclw_frac_alt_so[:,1],gfdl_hiram_tclw_frac_alt_so[:,0], '--g', label='CMIP5-GFDL-HIRAM-AMIP', linewidth=1)
ax1.plot(cam5_tclw_frac_alt_so[:,1],cam5_tclw_frac_alt_so[:,0], '--m', label='CMIP5-CESM1-CAM5-AMIP', linewidth=1)

ax1.plot(gfdl4_tclw_frac_alt_so[:,1],gfdl4_tclw_frac_alt_so[:,0], '-g', label='CMIP6-GFDL-AM4-AMIP', linewidth=1)
ax1.plot(cam6_tclw_frac_alt_so[:,1],cam6_tclw_frac_alt_so[:,0], '-m', label='CMIP6-CESM2-CAM6-AMIP', linewidth=1)

ax1.fill_betweenx(alt_lw_so, cmip5_alt_min, cmip5_alt_max, facecolor='black', alpha=0.2,
                label='CMIP5 Model Range')

ax1.fill_betweenx(alt_lw_so, cmip6_alt_min, cmip6_alt_max, facecolor='yellow', alpha=0.4,
                label='CMIP6 Model Range')



ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax1.set_ylabel('Altitude (km)')
ax1.set_xlabel('Cloud Liquid Water Fraction')

#ax1.set_title ('Global Cloud Fraction vs Latitude')
ax1.set_xlim(0, 0.6)
ax1.grid(True)
plt.savefig("tclw_frac_alt_so.pdf", format="pdf", bbox_inches='tight')
plt.show()
"""

#-----------TCLW_G Frac Temp-----------#
gfdl4_tclw_frac_temp_g = gfdl4_tclw_frac_temp_g[0:20]
temp_g = gfdl4_tclw_frac_temp_g[:,0]
gfdl4_tclw_frac_temp = gfdl4_tclw_frac_temp_g[:,1]

f = interpolate.interp1d(cccm_tclw_frac_temp_g[:,0], cccm_tclw_frac_temp_g[:,1], fill_value="extrapolate", kind = 'cubic')
cccm_tclw_frac_temp = f(temp_g)
cccm_tclw_frac_temp_g = np.vstack((temp_g, cccm_tclw_frac_temp)).T

f = interpolate.interp1d(calipso_tclw_frac_temp_g[:,0], calipso_tclw_frac_temp_g[:,1], fill_value="extrapolate", kind = 'cubic')
calipso_tclw_frac_temp = f(temp_g)
calipso_tclw_frac_temp_g = np.vstack((temp_g, calipso_tclw_frac_temp)).T

f = interpolate.interp1d(ecmwf_tclw_frac_temp_g[:,0], ecmwf_tclw_frac_temp_g[:,1], fill_value="extrapolate", kind = 'cubic')
ecmwf_tclw_frac_temp = f(temp_g)
ecmwf_tclw_frac_temp_g = np.vstack((temp_g, ecmwf_tclw_frac_temp)).T

f = interpolate.interp1d(gfdl_hiram_tclw_frac_temp_g[:,0], gfdl_hiram_tclw_frac_temp_g[:,1], fill_value="extrapolate")
gfdl_hiram_tclw_frac_temp = f(temp_g)
gfdl_hiram_tclw_frac_temp_g = np.vstack((temp_g, gfdl_hiram_tclw_frac_temp)).T

f = interpolate.interp1d(cam5_tclw_frac_temp_g[:,0], cam5_tclw_frac_temp_g[:,1], fill_value="extrapolate")
cam5_tclw_frac_temp = f(temp_g)
cam5_tclw_frac_temp_g = np.vstack((temp_g, cam5_tclw_frac_temp)).T

f = interpolate.interp1d(cam6_tclw_frac_temp_g[:,0], cam6_tclw_frac_temp_g[:,1], fill_value="extrapolate")
cam6_tclw_frac_temp = f(temp_g)
cam6_tclw_frac_temp_g = np.vstack((temp_g, cam6_tclw_frac_temp)).T

f = interpolate.interp1d(mri_cgcm_tclw_frac_temp_g[:,0], mri_cgcm_tclw_frac_temp_g[:,1], fill_value="extrapolate")
mri_cgcm_tclw_frac_temp = f(temp_g)
mri_cgcm_tclw_frac_temp_g = np.vstack((temp_g, mri_cgcm_tclw_frac_temp)).T

f = interpolate.interp1d(mri_tclw_frac_temp_g[:,0], mri_tclw_frac_temp_g[:,1], fill_value="extrapolate")
mri_tclw_frac_temp = f(temp_g)
mri_tclw_frac_temp_g = np.vstack((temp_g, mri_tclw_frac_temp)).T

f = interpolate.interp1d(miroc5_tclw_frac_temp_g[:,0], miroc5_tclw_frac_temp_g[:,1], fill_value="extrapolate")
miroc5_tclw_frac_temp = f(temp_g)
miroc5_tclw_frac_temp_g = np.vstack((temp_g, miroc5_tclw_frac_temp)).T

f = interpolate.interp1d(miroc6_tclw_frac_temp_g[:,0], miroc6_tclw_frac_temp_g[:,1], fill_value="extrapolate")
miroc6_tclw_frac_temp = f(temp_g)
miroc6_tclw_frac_temp_g = np.vstack((temp_g, miroc6_tclw_frac_temp)).T

f = interpolate.interp1d(ipsl5_tclw_frac_temp_g[:,0], ipsl5_tclw_frac_temp_g[:,1], fill_value="extrapolate")
ipsl5_tclw_frac_temp = f(temp_g)
ipsl5_tclw_frac_temp_g = np.vstack((temp_g, ipsl5_tclw_frac_temp)).T

f = interpolate.interp1d(ipsl6_tclw_frac_temp_g[:,0], ipsl6_tclw_frac_temp_g[:,1], fill_value="extrapolate")
ipsl6_tclw_frac_temp = f(temp_g)
ipsl6_tclw_frac_temp_g = np.vstack((temp_g, ipsl6_tclw_frac_temp)).T

f = interpolate.interp1d(giss5_tclw_frac_temp_g[:,0], giss5_tclw_frac_temp_g[:,1], fill_value="extrapolate")
giss5_tclw_frac_temp = f(temp_g)
giss5_tclw_frac_temp_g = np.vstack((temp_g, giss5_tclw_frac_temp)).T

f = interpolate.interp1d(giss6_tclw_frac_temp_g[:,0], giss6_tclw_frac_temp_g[:,1], fill_value="extrapolate")
giss6_tclw_frac_temp = f(temp_g)
giss6_tclw_frac_temp_g = np.vstack((temp_g, giss6_tclw_frac_temp)).T


cmip5_alt_max = np.maximum.reduce([ecmwf_tclw_frac_temp, gfdl_hiram_tclw_frac_temp, cam5_tclw_frac_temp, mri_cgcm_tclw_frac_temp, miroc5_tclw_frac_temp, ipsl5_tclw_frac_temp, giss5_tclw_frac_temp])
cmip5_alt_min = np.minimum.reduce([ecmwf_tclw_frac_temp, gfdl_hiram_tclw_frac_temp, cam5_tclw_frac_temp, mri_cgcm_tclw_frac_temp, miroc5_tclw_frac_temp, ipsl5_tclw_frac_temp, giss5_tclw_frac_temp])

cmip6_alt_max = np.maximum.reduce([ecmwf_tclw_frac_temp, gfdl4_tclw_frac_temp, cam6_tclw_frac_temp, mri_tclw_frac_temp, miroc6_tclw_frac_temp, ipsl6_tclw_frac_temp, giss6_tclw_frac_temp])
cmip6_alt_min = np.minimum.reduce([ecmwf_tclw_frac_temp, gfdl4_tclw_frac_temp, cam6_tclw_frac_temp, mri_tclw_frac_temp, miroc6_tclw_frac_temp, ipsl6_tclw_frac_temp, giss6_tclw_frac_temp])

"""
fig, ax1 = plt.subplots()


ax1.plot(calipso_tclw_frac_temp_g[:,0],calipso_tclw_frac_temp_g[:,1], '-b', label='CALIPSO-GOCCP', linewidth=1)
ax1.plot(cccm_tclw_frac_temp_g[:,0],cccm_tclw_frac_temp_g[:,1], '--b', label='CCCM - Satellite', linewidth=1)

ax1.plot(gfdl_hiram_tclw_frac_temp_g[:,0],gfdl_hiram_tclw_frac_temp_g[:,1], '--g', label='CMIP5-GFDL-HIRAM-AMIP', linewidth=1)
ax1.plot(cam5_tclw_frac_temp_g[:,0],cam5_tclw_frac_temp_g[:,1], '--m', label='CMIP5-CESM1-CAM5-AMIP', linewidth=1)


ax1.plot(gfdl4_tclw_frac_temp_g[:,0],gfdl4_tclw_frac_temp_g[:,1], '-g', label='CMIP6-GFDL-AM4-AMIP', linewidth=1)
ax1.plot(cam6_tclw_frac_temp_g[5:,0],cam6_tclw_frac_temp_g[5:,1], '-m', label='CMIP6-CESM2-CAM6-AMIP', linewidth=1)


ax1.fill_between(temp_g, cmip5_alt_min, cmip5_alt_max, facecolor='black', alpha=0.4,
                label='CMIP5 Model Range')

ax1.fill_between(temp_g, cmip6_alt_min, cmip6_alt_max, facecolor='yellow', alpha=0.4,
                label='CMIP6 Model Range')


ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));


ax1.set_ylabel('Liquid Fraction')
ax1.set_xlabel('Temperature (K)')

#ax1.set_title ('Global Cloud Fraction vs Latitude')

ax1.grid(True)
plt.savefig("tclw_frac_temp_g.pdf", format="pdf", bbox_inches='tight')
plt.show()
"""

#-----------TCLW_SO Frac Temp-----------#
gfdl4_tclw_frac_temp_so = gfdl4_tclw_frac_temp_so[0:20]
temp_so = gfdl4_tclw_frac_temp_so[:,0]
gfdl4_tclw_frac_temp_s = gfdl4_tclw_frac_temp_so[:,1]

f = interpolate.interp1d(cccm_tclw_frac_temp_so[:,0], cccm_tclw_frac_temp_so[:,1], fill_value="extrapolate", kind = 'cubic')
cccm_tclw_frac_temp_s = f(temp_so)
cccm_tclw_frac_temp_so = np.vstack((temp_so, cccm_tclw_frac_temp_s)).T

f = interpolate.interp1d(calipso_tclw_frac_temp_so[:,0], calipso_tclw_frac_temp_so[:,1], fill_value="extrapolate")
calipso_tclw_frac_temp_s = f(temp_so)
calipso_tclw_frac_temp_so = np.vstack((temp_so, calipso_tclw_frac_temp_s)).T

f = interpolate.interp1d(ecmwf_tclw_frac_temp_so[:,0], ecmwf_tclw_frac_temp_so[:,1], fill_value="extrapolate", kind = 'cubic')
ecmwf_tclw_frac_temp_s = f(temp_so)
ecmwf_tclw_frac_temp_so = np.vstack((temp_so, ecmwf_tclw_frac_temp_s)).T

f = interpolate.interp1d(gfdl_hiram_tclw_frac_temp_so[:,0], gfdl_hiram_tclw_frac_temp_so[:,1], fill_value="extrapolate")
gfdl_hiram_tclw_frac_temp_s = f(temp_so)
gfdl_hiram_tclw_frac_temp_so = np.vstack((temp_so, gfdl_hiram_tclw_frac_temp_s)).T

f = interpolate.interp1d(cam5_tclw_frac_temp_so[:,0], cam5_tclw_frac_temp_so[:,1], fill_value="extrapolate")
cam5_tclw_frac_temp_s = f(temp_so)
cam5_tclw_frac_temp_so = np.vstack((temp_so, cam5_tclw_frac_temp_s)).T

f = interpolate.interp1d(cam6_tclw_frac_temp_so[:,0], cam6_tclw_frac_temp_so[:,1], fill_value="extrapolate")
cam6_tclw_frac_temp_s = f(temp_so)
cam6_tclw_frac_temp_so = np.vstack((temp_so, cam6_tclw_frac_temp_s)).T

f = interpolate.interp1d(mri_cgcm_tclw_frac_temp_so[:,0], mri_cgcm_tclw_frac_temp_so[:,1], fill_value="extrapolate")
mri_cgcm_tclw_frac_temp_s = f(temp_so)
mri_cgcm_tclw_frac_temp_so = np.vstack((temp_so, mri_cgcm_tclw_frac_temp_s)).T

f = interpolate.interp1d(mri_tclw_frac_temp_so[:,0], mri_tclw_frac_temp_so[:,1], fill_value="extrapolate")
mri_tclw_frac_temp_s = f(temp_so)
mri_tclw_frac_temp_so = np.vstack((temp_so, mri_tclw_frac_temp_s)).T

f = interpolate.interp1d(miroc5_tclw_frac_temp_so[:,0], miroc5_tclw_frac_temp_so[:,1], fill_value="extrapolate")
miroc5_tclw_frac_temp_s = f(temp_so)
miroc5_tclw_frac_temp_so = np.vstack((temp_so, miroc5_tclw_frac_temp_s)).T

f = interpolate.interp1d(miroc6_tclw_frac_temp_so[:,0], miroc6_tclw_frac_temp_so[:,1], fill_value="extrapolate")
miroc6_tclw_frac_temp_s = f(temp_so)
miroc6_tclw_frac_temp_so = np.vstack((temp_so, miroc6_tclw_frac_temp_s)).T

f = interpolate.interp1d(ipsl5_tclw_frac_temp_so[:,0], ipsl5_tclw_frac_temp_so[:,1], fill_value="extrapolate")
ipsl5_tclw_frac_temp_s = f(temp_so)
ipsl5_tclw_frac_temp_so = np.vstack((temp_so, ipsl5_tclw_frac_temp_s)).T

f = interpolate.interp1d(ipsl6_tclw_frac_temp_so[:,0], ipsl6_tclw_frac_temp_so[:,1], fill_value="extrapolate")
ipsl6_tclw_frac_temp_s = f(temp_so)
ipsl6_tclw_frac_temp_so = np.vstack((temp_so, ipsl6_tclw_frac_temp_s)).T

f = interpolate.interp1d(giss5_tclw_frac_temp_so[:,0], giss5_tclw_frac_temp_so[:,1], fill_value="extrapolate")
giss5_tclw_frac_temp_s = f(temp_so)
giss5_tclw_frac_temp_so = np.vstack((temp_so, giss5_tclw_frac_temp_s)).T

f = interpolate.interp1d(giss6_tclw_frac_temp_so[:,0], giss6_tclw_frac_temp_so[:,1], fill_value="extrapolate")
giss6_tclw_frac_temp_s = f(temp_so)
giss6_tclw_frac_temp_so = np.vstack((temp_so, giss6_tclw_frac_temp_s)).T


cmip5_alt_max = np.maximum.reduce([ecmwf_tclw_frac_temp_s, gfdl_hiram_tclw_frac_temp_s, cam5_tclw_frac_temp_s, mri_cgcm_tclw_frac_temp_s, miroc5_tclw_frac_temp_s, ipsl5_tclw_frac_temp_s, giss5_tclw_frac_temp_s])
cmip5_alt_min = np.minimum.reduce([ecmwf_tclw_frac_temp_s, gfdl_hiram_tclw_frac_temp_s, cam5_tclw_frac_temp_s, mri_cgcm_tclw_frac_temp_s, miroc5_tclw_frac_temp_s, ipsl5_tclw_frac_temp_s, giss5_tclw_frac_temp_s])

cmip6_alt_max = np.maximum.reduce([ecmwf_tclw_frac_temp_s, gfdl4_tclw_frac_temp_s, cam6_tclw_frac_temp_s, mri_tclw_frac_temp_s, miroc6_tclw_frac_temp_s, ipsl6_tclw_frac_temp_s, giss6_tclw_frac_temp_s])
cmip6_alt_min = np.minimum.reduce([ecmwf_tclw_frac_temp_s, gfdl4_tclw_frac_temp_s, cam6_tclw_frac_temp_s, mri_tclw_frac_temp_s, miroc6_tclw_frac_temp_s, ipsl6_tclw_frac_temp_s, giss6_tclw_frac_temp_s])

"""
fig, ax1 = plt.subplots()


ax1.plot(calipso_tclw_frac_temp_so[:,0],calipso_tclw_frac_temp_so[:,1], '-b', label='CALIPSO-GOCCP', linewidth=1)
ax1.plot(cccm_tclw_frac_temp_so[:,0],cccm_tclw_frac_temp_so[:,1], '--b', label='CCCM - Satellite', linewidth=1)

ax1.plot(gfdl_hiram_tclw_frac_temp_so[:,0],gfdl_hiram_tclw_frac_temp_so[:,1], '--g', label='CMIP5-GFDL-HIRAM-AMIP', linewidth=1)
ax1.plot(cam5_tclw_frac_temp_so[:,0],cam5_tclw_frac_temp_so[:,1], '--m', label='CMIP5-CESM1-CAM5-AMIP', linewidth=1)


ax1.plot(gfdl4_tclw_frac_temp_so[:,0],gfdl4_tclw_frac_temp_so[:,1], '-g', label='CMIP6-GFDL-AM4-AMIP', linewidth=1)
ax1.plot(cam6_tclw_frac_temp_so[5:,0],cam6_tclw_frac_temp_so[5:,1], '-m', label='CMIP6-CESM2-CAM6-AMIP', linewidth=1)


ax1.fill_between(temp_so, cmip5_alt_min, cmip5_alt_max, facecolor='black', alpha=0.4,
                label='CMIP5 Model Range')

ax1.fill_between(temp_so, cmip6_alt_min, cmip6_alt_max, facecolor='yellow', alpha=0.4,
                label='CMIP6 Model Range')


ax1.axvline(x=273, label = '273K', color = 'black', linestyle='--')
ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));


ax1.set_ylabel('Liquid Fraction')
ax1.set_xlabel('Temperature (K)')

#ax1.set_title ('Global Cloud Fraction vs Latitude')

ax1.grid(True)
plt.savefig("tclw_frac_temp_so.pdf", format="pdf", bbox_inches='tight')
plt.show()
"""
