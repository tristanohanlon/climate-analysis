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
        
grid_data_sets = [
        grid_DataSet('Cloud free area percent coverage (CALIPSO-CloudSat)' ),
        grid_DataSet('Liquid water content profile used' ),
        grid_DataSet('Ice water content profile used' ),
]        

# The directory where your HDF files are stored
os.chdir('E:/CCCM/test')  # Home PC



# Load every file in the directory
for filename in os.listdir(): 
    pprint( filename )
    # Load the file
    f = SD.SD(filename)
    raw_lon = f.select('Longitude of subsatellite point at surface at observation').get()
    raw_lat = 90 - f.select('Colatitude of subsatellite point at surface at observation').get()
    for data_set in grid_data_sets:
        sel = f.select( data_set.type_name )
        data = sel.get()
        if data_set.type_name == 'Cloud free area percent coverage (CALIPSO-CloudSat)':
            data = (100 - data) / 100
        else:
            fill_value = sel.attributes()['_FillValue']
            data[ data == fill_value ] = None #set fill values to nan
            data = np.nansum(data, axis = 1)
        for l_index, ( la, lo ) in enumerate( zip( raw_lat, raw_lon ) ):
            if la <= constants.min_lat or la >= constants.max_lat:
                continue
            
            lat_bin = int( ( la - constants.min_lat ) / constants.lat_division)
            lon_bin = int( lo - 1 )
            
            #print( lat_bin, alt_bin )
            val = data[l_index]
            if val == None:
                continue
            data_set.data[ lat_bin, lon_bin ] += val
            data_set.data_counts[ lat_bin, lon_bin ] += 1
            
for data_set in grid_data_sets:
    data_set.data /= data_set.data_counts
#    pprint( data_set.data.shape )
  
clt_lat_lon = data_sets[0].data
clwvi_lat_lon = data_sets[1].data
clivi_lat_lon = data_sets[2].data

#---create liquid and ice fractions---#

lw_clwvi_lat_lon = (clwvi_lat_lon / (clwvi_lat_lon + clivi_lat_lon)) * clt_lat_lon
iw_clivi_lat_lon = (clivi_lat_lon / (clwvi_lat_lon + clivi_lat_lon)) * clt_lat_lon


clt = np.nanmean( clt_lat_lon, axis = -1 )
clwvi = np.nanmean( lw_clwvi_lat_lon, axis = -1 )
clivi = np.nanmean( iw_clivi_lat_lon, axis = -1 )

#----------------------------#

os.chdir( location + '/climate-analysis/reduced_data' )

f = h5py.File('Jun_2006_Apr_2011_CCCM.h5','a')
f.create_dataset('clt', data=clt)
f.create_dataset('clwvi', data=clwvi)
f.create_dataset('clivi', data=clivi)
f.create_dataset('clt_lat_lon', data=clt_lat_lon)
f.create_dataset('clwvi_lat_lon', data=lw_clwvi_lat_lon)
f.create_dataset('clivi_lat_lon', data=iw_clivi_lat_lon)
f.close()

































