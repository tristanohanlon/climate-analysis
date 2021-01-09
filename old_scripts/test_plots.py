# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 09:47:57 2019

@author: Tristan
"""
import h5py
import os
import matplotlib.pyplot as plt
import constants
import numpy as np
import cartopy.crs as ccrs
from netCDF4 import Dataset
from netCDF4 import date2index
import math
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline
import datetime
from pprint import pprint
import constants

def extract_data_over_time( type, f ):
        data = np.array( f.variables[type][:] )
        return data

model = 'CMIP6-CESM2-CAM6' 
location = constants.home # home, uni, hdd or laptop
os.chdir( location + 'Data/' + model )

with Dataset( 'ncl_a1_SRF' + constants.model_dict[ model ], 'r') as f:
    ncl_a1_SRF = extract_data_over_time('ncl_a1_SRF', f )
    lon = f.variables['lon'][:]
    lat = f.variables['lat'][:]
ncl_a1_SRF = np.nanmean(ncl_a1_SRF, axis = 0)


with Dataset( 'ncl_a2_SRF' + constants.model_dict[ model ], 'r') as f:
    ncl_a2_SRF = extract_data_over_time('ncl_a2_SRF', f )
ncl_a2_SRF = np.nanmean(ncl_a2_SRF, axis = 0)

with Dataset( 'ncl_a3_SRF' + constants.model_dict[ model ], 'r') as f:
    ncl_a3_SRF = extract_data_over_time('ncl_a3_SRF', f )
ncl_a3_SRF = np.nanmean(ncl_a3_SRF, axis = 0)


with Dataset( 'num_a4_SRF' + constants.model_dict[ model ], 'r') as f:
    num_a4_SRF = extract_data_over_time('num_a4_SRF', f )
num_a4_SRF = np.nanmean(num_a4_SRF, axis = 0)

num = ncl_a1_SRF + ncl_a2_SRF+ ncl_a3_SRF

ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
ax.coastlines()
p = ax.contourf(lon, lat, num, transform=ccrs.PlateCarree(), cmap='coolwarm')
cbar = plt.colorbar(p, orientation='horizontal')
ax.set_title ('Total Normalised Aerosol Cover (Dust, Organic, S04)')
#cbar.set_clim(-45, 45)

cbar.set_label('Total Normalised Aerosol Cover (Dust, Organic, S04)')
plt.show()
