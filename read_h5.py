# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 14:44:25 2019

@author: Tristan

    CMIP5-CESM1-CAM5
    CMIP5-GFDL-HIRAM-C360
    CMIP5-GISS-E2R
    CMIP5-IPSL-CM5A-LR
    CMIP5-MIROC5
    CMIP5-MRI-CGCM3
    
    CMIP6-CESM2-CAM6
    CMIP6-GFDL-AM4
    CMIP6-GISS-E21G
    CMIP6-IPSL-CM6A-LR
    CMIP6-MIROC6
    CMIP6-MRI-ESM2
    
    CCCM
    ECMWF
"""

#read hdf5
import h5py
import os
import matplotlib.pyplot as plt
import constants
import numpy as np

#--- Set Location, date period and model ---#

model = 'CCCM'
location = constants.home + '/climate-analysis/reduced_data'
date_cmip5 = 'Jan_2001_Dec_2005'
date_cmip6 = 'Jan_2006_Dec_2010'
os.chdir( location )

if model == 'CCCM' or model == 'ECMWF':
    h5f = h5py.File( model + '.h5', 'r')
    cl_alt_lat = h5f['cl_alt_lat'][:]
    cl_g = h5f['cl_g'][:]
    cl_so = h5f['cl_so'][:]
    cli_alt_lat = h5f['cli_alt_lat'][:]
    cli_g = h5f['cli_g'][:]
    cli_so = h5f['cli_so'][:]
    clw_alt_lat = h5f['clw_alt_lat'][:]
    clw_g = h5f['clw_g'][:]
    clw_so = h5f['clw_so'][:]
    ta_alt_lat = h5f['ta_alt_lat'][:]
    ta_liq_alt_lat = h5f['ta_liq_alt_lat'][:]
    ta_liq_g = h5f['ta_liq_g'][:]
    ta_liq_so = h5f['ta_liq_so'][:]
    lw_frac_alt_lat = h5f['lw_frac_alt_lat'][:]

else:
    h5f = h5py.File( date_cmip6 + '_' + model + '.h5', 'r')
    cl_alt_lat = h5f['cl_alt_lat'][:]
    cl_g = h5f['cl_g'][:]
    cl_so = h5f['cl_so'][:]
    cli_alt_lat = h5f['cli_alt_lat'][:]
    cli_g = h5f['cli_g'][:]
    cli_so = h5f['cli_so'][:]
    clt = h5f['clt'][:]
    clw_alt_lat = h5f['clw_alt_lat'][:]
    clw_g = h5f['clw_g'][:]
    clw_so = h5f['clw_so'][:]
    ta_alt_lat = h5f['ta_alt_lat'][:]
    ta_liq_g = h5f['ta_liq_g'][:]
    ta_liq_so = h5f['ta_liq_so'][:]

for index,key in enumerate(h5f.keys()):
    print (index, key)
    
#--- sample plots for confirmation ---#

fig, ax = plt.subplots()
ax.plot( cl_so, constants.alt )
plt.grid(True)
plt.show()

fig, ax = plt.subplots()
ax.plot( ta_liq_g, clw_g )
plt.grid(True)
plt.show()

fig, ax = plt.subplots()
cont = ax.contourf( constants.lat, constants.liq_alt, clw_alt_lat )
temp = ax.contour( constants.lat, constants.liq_alt, (ta_liq_alt_lat - 273.15), colors='white')
temp.collections[5].set_linewidth(3)
temp.collections[5].set_color('white')
ax.clabel(temp, inline=1, fontsize=10)

cbar = fig.colorbar(cont)
plt.show()
