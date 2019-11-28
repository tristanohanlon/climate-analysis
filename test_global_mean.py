"""
Created on Sun Oct  6 14:44:25 2019

@author: Tristan O'Hanlon - University of Auckland

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
    CERES
    CALIPSO
"""

import h5py
import os
import matplotlib.pyplot as plt
import constants
import numpy as np
import cartopy.crs as ccrs
import datetime
from netCDF4 import Dataset
from netCDF4 import date2index

#--- Set Location, date period and model ---#

# specify model from the list above
# specify location: home, uni, hdd, laptop
model = 'CMIP6-GFDL-AM4' # see comments above for available models
location = constants.hdd # home, uni, hdd or laptop
os.chdir( location + 'Data/' + model )
variable = 'cl'


start = datetime.datetime( 2006, 1, 1 )
end = datetime.datetime( 2010, 12, 1 )


with Dataset( variable + constants.model_dict[ model ], 'r') as f:
        raw_lat = constants.extract_data( f, 'lat')
        raw_lon = constants.extract_data( f, 'lon')
        raw_lev = constants.extract_data( f, 'lev')

        data = np.nanmean( constants.extract_data_over_time(variable, f, start, end ), axis = 0 ) # average over time
        data = constants.global3DMean(data, raw_lat)
        print(data)


fig, ax = plt.subplots()
ax.plot( data, raw_lev )
ax.set_ylabel('Altitude (km)')
plt.grid(True)
plt.show()





