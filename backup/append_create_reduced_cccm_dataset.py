# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon - University of Auckland, Samual Crookes & Jonathan Rogers


"""
from pprint import pprint
import time
import numpy as np
import os
from pyhdf import SD
import h5py
import matplotlib.pyplot as plt
from scipy import interpolate
import constants
from scipy import stats
import cartopy.crs as ccrs

###############################################################################
start = time.time()

#set location

location = constants.home
#lat = constants.lat


class grid_DataSet:
    def __init__( self, type_name ):
        self.type_name = type_name
        self.data = np.zeros(( constants.lat.size, constants.lon.size ))
        self.data_counts = np.zeros( (constants.lat.size, constants.lon.size ))


altitude_types = [ 'Layer center height profile (clouds and aerosol)' ]

class DataSet:
    def __init__( self, type_name ):
        self.type_name = type_name
        self.data = np.zeros(( constants.lat.size, constants.alt.size ))
        self.data_counts = np.zeros( (constants.lat.size, constants.alt.size ))
        

data_sets = [
        # DataSet('Cloud fraction profile' ),
        DataSet('Liquid water content profile used' ), # kg/m3
        DataSet('Ice water content profile used' ), # kg/m3
        # DataSet('Temperature profile' )        
]        

grid_data_sets = [
        grid_DataSet('Cloud free area percent coverage (CALIPSO-CloudSat)' ),
        grid_DataSet('Liquid water content profile used' ),
        grid_DataSet('Ice water content profile used' ),
]        

# The directory where your HDF files are stored
os.chdir('E:/CCCM/test_1')  # Home PC


# Load every file in the directory
for filename in os.listdir(): 
    # Load the file
    if os.path.isdir( filename ):
        continue
    pprint( filename )
    # Load the file
    f = SD.SD(filename)
    raw_lon = f.select('Longitude of subsatellite point at surface at observation').get()
    raw_lat = 90 - f.select('Colatitude of subsatellite point at surface at observation').get()
    raw_alt_phase = f.select('Irradiance layer center height profile').get() / 1000
    # altitudes[0] /= 1000
    # altitudes[2] /= 1000

    for data_set in grid_data_sets:
        sel = f.select( data_set.type_name )
        data = sel.get()
        if data_set.type_name == 'Cloud free area percent coverage (CALIPSO-CloudSat)':
            data = (100 - data) / 100
        else:
            fill_value = sel.attributes()['_FillValue']
            data[ data == fill_value ] = None #set fill values to nan
            data = np.nanmean(data, axis = 1)


        data = stats.binned_statistic_2d(raw_lat, raw_lon, data, np.nanmean, bins=[constants.lat.size,constants.lon.size], range=[[constants.min_lat, constants.max_lat], [constants.min_lon, constants.max_lon]])
        data_count = stats.binned_statistic_2d(raw_lat, raw_lon, data, 'count', bins=[constants.lat.size,constants.lon.size], range=[[constants.min_lat, constants.max_lat], [constants.min_lon, constants.max_lon]])

        data_set.data = data_set.data + data.statistic
        data_set.data_counts = data_set.data_counts + data_count.statistic


    for data_set in data_sets:
        sel = f.select( data_set.type_name )
        data = np.transpose(sel.get()).tolist()

        fill_value = sel.attributes()['_FillValue']
        data[ data == fill_value ] = None #set fill values to nan

        data = stats.binned_statistic_2d(raw_lat, raw_alt_phase, data, np.nanmean, bins=[constants.lat.size,constants.alt.size], range=[[constants.min_lat, constants.max_lat], [constants.min_alt, constants.max_alt]])
        data_count = stats.binned_statistic_2d(raw_lat, raw_alt_phase, data, 'count', bins=[constants.lat.size,constants.alt.size], range=[[constants.min_lat, constants.max_lat], [constants.min_alt, constants.max_alt]])

        data_set.data = data_set.data + data.statistic
        data_set.data_counts = data_set.data_counts + data_count.statistic


for data_set in grid_data_sets:
    data_set.data /= data_set.data_counts
#    pprint( data_set.data.shape )
  
clt_lat_lon = grid_data_sets[0].data
clwvi_lat_lon = grid_data_sets[1].data
clivi_lat_lon = grid_data_sets[2].data

# for data_set in data_sets:
#     data_set.data /= data_set.data_counts
# #    pprint( data_set.data.shape )
  
# cl_alt_lat = data_sets[0].data
# clw_alt_lat = data_sets[1].data
# cli_alt_lat = data_sets[2].data
# full_ta_alt_lat = data_sets[3].data

#---create liquid and ice fractions---#


clt = np.nanmean( clt_lat_lon, axis = -1 )
clwvi = np.nanmean( clwvi_lat_lon, axis = -1 )
clivi = np.nanmean( clivi_lat_lon, axis = -1 )


fig, ax = plt.subplots()
ax.plot( constants.lat[constants.lat_confine_1:constants.lat_confine_2], clt[constants.lat_confine_1:constants.lat_confine_2] )
ax.set_ylabel('Cloud Fraction')
ax.set_xlabel('Latitude')
ax.set_title ('Global Cloud Fraction vs Latitude')
plt.grid(True)
plt.show()

ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
ax.coastlines()
p = ax.contourf(constants.lon, constants.lat, clt_lat_lon, transform=ccrs.PlateCarree(), cmap='coolwarm')
cbar = plt.colorbar(p, orientation='horizontal')
cbar.set_label('Cloud Fraction')
plt.show()

end = time.time()
print('Create dataset took:', end - start, 's')
#----------------------------#

os.chdir( location + '/climate-analysis/reduced_data' )

f = h5py.File('test_CCCM.h5','a')

f.create_dataset('clt', data=clt)
f.create_dataset('clwvi', data=clwvi)
f.create_dataset('clivi', data=clivi)
f.create_dataset('clt_lat_lon', data=clt_lat_lon)
f.create_dataset('clwvi_lat_lon', data=lw_clwvi_lat_lon)
f.create_dataset('clivi_lat_lon', data=iw_clivi_lat_lon)
f.close()

































