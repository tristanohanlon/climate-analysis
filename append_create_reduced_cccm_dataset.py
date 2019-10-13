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


class DataSet:
    def __init__( self, type_name ):
        self.type_name = type_name
        self.data = np.zeros(( constants.lat.size, constants.lon.size ))
        self.data_counts = np.zeros( (constants.lat.size, constants.lon.size ))
        
data_sets = [
        DataSet('Cloud fraction profile' ),
        DataSet('Liquid water content profile used' ),
        DataSet('Ice water content profile used' ),
]        

# The directory where your HDF files are stored
os.chdir('E:/CCCM/2007')  # Home PC



# Load every file in the directory
for filename in os.listdir(): 
    pprint( filename )
    # Load the file
    f = SD.SD(filename)
    raw_lon = f.select('Longitude of subsatellite point at surface at observation').get()
    raw_lat = 90 - f.select('Colatitude of subsatellite point at surface at observation').get()
    for data_set in data_sets:
        sel = f.select( data_set.type_name )
        data = sel.get()
        data = np.nansum(data, axis = 1)
        fill_value = sel.attributes()['_FillValue']
        for l_index, l in enumerate( raw_lat ):
            if l <= constants.min_lat or l >= constants.max_lat:
                continue
            lat_bin = int( ( l - constants.min_lat ) / constants.lat_division)
            for a_index, a in enumerate( raw_lon ):
                lon_bin = int( ( a - constants.min_lon ) / constants.lon_division )
                #print( lat_bin, alt_bin )
                val = data[l_index, a_index]
                if val == fill_value:
                    continue
                data_set.data[ lat_bin, lon_bin ] += val
                data_set.data_counts[ lat_bin, lon_bin ] += 1
            
for data_set in data_sets:
    data_set.data /= data_set.data_counts
#    pprint( data_set.data.shape )
  
clt_lat_lon = data_sets[0].data
clwvi_lat_lon = data_sets[1].data
clivi_lat_lon = data_sets[2].data

#---create liquid and ice fractions---#

lw_clwvi_lat_lon = (clwvi_lat_lon / (clwvi_lat_lon + clivi_lat_lon)) * clt_lat_lon
cli_alt_lat = (clivi_lat_lon / (clwvi_lat_lon + clivi_lat_lon)) * clt_lat_lon
clwvi_lat_lon = lw_clwvi_lat_lon

#----------------------------#

#os.chdir( location + '/climate-analysis/reduced_data' )
#
#save_filename = 'Jun_2006_Apr_2011_CCCM.h5'
#
#with h5py.File(save_filename, 'w') as p:
#            
#    p.create_dataset('cl_g', data=cl_g) # global layer total cloud fraction corresponding to alt
#    p.create_dataset('clw_g', data=clw_g) # global layer cloud liquid water fraction corresponding to liq_alt
#    p.create_dataset('cli_g', data=cli_g) # global layer cloud ice water fraction corresponding to alt
#    
#    p.create_dataset('cl_so', data=cl_so) # southern ocean layer total cloud fraction corresponding to alt
#    p.create_dataset('clw_so', data=clw_so) # southern ocean layer cloud liquid water fraction corresponding to liq_alt
#    p.create_dataset('cli_so', data=cli_so) # southern ocean layer cloud ice water fraction corresponding to alt
#    
#    p.create_dataset('clw_t_g', data=clw_t_g) # global layer cloud liquid water fraction corresponding to ta_g
#    p.create_dataset('clw_t_so', data=clw_t_so) # global layer cloud liquid water fraction corresponding to ta_so
#    
#    p.create_dataset('ta_alt_lat', data= np.transpose( ta_alt_lat ) ) # temperature corresponding to liq_alt and lat
#    p.create_dataset('cl_alt_lat', data= np.transpose( cl_alt_lat ) ) # total cloud fraction corresponding to alt and lat
#    p.create_dataset('clw_alt_lat', data= np.transpose( clw_alt_lat ) ) # cloud liquid water fraction corresponding to liq_alt and lat
#    p.create_dataset('cli_alt_lat', data= np.transpose( cli_alt_lat ) ) # cloud ice water fraction corresponding to alt and lat
#
#    p.create_dataset('full_ta_alt_lat', data= np.transpose( full_ta_alt_lat ) ) # temperature corresponding to alt and lat
#    p.create_dataset('full_clw_alt_lat', data= np.transpose( full_clw_alt_lat ) ) # tcloud liquid water fraction corresponding to alt and lat
#
#    p.close()
#
#end = time.time()
#print('Create dataset took:', end - start, 's')



































