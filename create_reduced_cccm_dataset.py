# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon


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

altitude_types = [ 'Irradiance layer center height profile', 'Layer center height profile (clouds and aerosol)', 'Irradiance level height profile' ]

class DataSet:
    def __init__( self, type_name, altitude_type ):
        self.type_name = type_name
        self.altitude_type = altitude_type
        self.data = np.zeros(( constants.lat.size, constants.alt.size ))
        self.data_counts = np.zeros( (constants.lat.size, constants.alt.size ))
        
data_sets = [
        DataSet('Cloud fraction profile', 1 ),
        DataSet('Liquid water content profile used', 0 ),
        DataSet('Ice water content profile used', 0 ),
        DataSet('Temperature profile', 2 )        
]        

# The directory where your HDF files are stored
os.chdir('E:/CCCM/test')  # Home PC

def create_southern_ocean_data( lat, global_data ):
    
    start_index = 0
    end_index = 0
    for index, l in enumerate( lat ):
        if l < -70:
            start_index = index + 1
        if l <= -50:
            end_index = index
    so = global_data[start_index:end_index]
    so = np.nanmean( so, axis = 0 ) # average over latitude
    return so


# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    
    raw_lat = f.select('Colatitude of subsatellite point at surface at observation').get() - 90
    
    altitudes = []
    for altitude_type in altitude_types:
        altitudes.append( f.select(altitude_type).get() )
    altitudes[0] /= 1000
    altitudes[2] /= 1000
    
    for data_set in data_sets:
        sel = f.select( data_set.type_name )
        data = sel.get()
        fill_value = sel.attributes()['_FillValue']
        for l_index, l in enumerate( raw_lat ):
            if l <= constants.min_lat or l >= constants.max_lat:
                continue
            lat_bin = int( ( l - constants.min_lat ) / constants.lat_division)
            for a_index, a in enumerate( altitudes[data_set.altitude_type] ):
                if a <= constants.min_alt or a >= constants.max_alt:
                    continue
                alt_bin = int( ( a - constants.min_alt ) / constants.alt_division )
                #print( lat_bin, alt_bin )
                val = data[l_index, a_index]
                if val == fill_value:
                    continue
                data_set.data[ lat_bin, alt_bin ] += val
                data_set.data_counts[ lat_bin, alt_bin ] += 1
            
for data_set in data_sets:
    data_set.data /= data_set.data_counts
#    pprint( data_set.data.shape )
  
cl_alt_lat = data_sets[0].data
clw_alt_lat = data_sets[1].data
cli_alt_lat = data_sets[2].data
ta_alt_lat = data_sets[3].data

#---create liquid and ice fractions---#

lw_frac_alt_lat = (clw_alt_lat / (clw_alt_lat + cli_alt_lat)) * cl_alt_lat
cli_alt_lat = (cli_alt_lat / (clw_alt_lat + cli_alt_lat)) * cl_alt_lat

#---create reduced altitude liquid fraction and temperatures---#

new_liq_data = np.zeros(( constants.lat.size, constants.liq_alt.size ))
liq_data_counts = np.zeros( (constants.lat.size, constants.liq_alt.size ))

for l_index, l in enumerate( constants.lat ):
    if l < constants.min_lat or l > constants.max_lat:
        continue
    lat_bin = int( ( l - constants.min_lat ) / constants.lat_division)
    for a_index, a in enumerate( constants.liq_alt ):
        if a < constants.min_liq_alt or a > constants.max_liq_alt:
            continue
        liq_alt_bin = int( ( a - constants.min_liq_alt ) / constants.liq_alt_division )
        val = lw_frac_alt_lat[l_index, a_index]
        if np.isnan(val):
            continue
        new_liq_data[ lat_bin, liq_alt_bin ] += val
        liq_data_counts[ lat_bin, liq_alt_bin ] += 1
            
clw_alt_lat = new_liq_data / liq_data_counts

interpolated = interpolate.interp2d( constants.alt, constants.lat, ta_alt_lat, kind = 'cubic', fill_value="extrapolate")
ta_liq_alt_lat = interpolated( constants.liq_alt, constants.lat )

#---create southern ocean and global datasets---#
    
cl_so = create_southern_ocean_data( constants.lat, cl_alt_lat)
clw_so = create_southern_ocean_data( constants.lat, clw_alt_lat)
cli_so = create_southern_ocean_data( constants.lat, cli_alt_lat)
ta_liq_so = create_southern_ocean_data( constants.lat, ta_liq_alt_lat)

cl_g = np.nanmean( cl_alt_lat, axis = 0 )
clw_g = np.nanmean( clw_alt_lat, axis = 0 )
cli_g = np.nanmean( cli_alt_lat, axis = 0 )
ta_liq_g = np.nanmean( ta_liq_alt_lat, axis = 0 )

#----------------------------#

os.chdir( location + '/climate-analysis/reduced_data' )

save_filename = 'CCCM.h5'

with h5py.File(save_filename, 'w') as p:
    
    p.create_dataset('ta_liq_g', data=ta_liq_g) # global layer temperature corresponding to liq_alt
    p.create_dataset('ta_liq_so', data=ta_liq_so) # southern ocean layer temperature corresponding to liq_alt
         
    p.create_dataset('cl_g', data=cl_g) # global layer total cloud fraction corresponding to alt
    p.create_dataset('clw_g', data=clw_g) # global layer cloud liquid water fraction corresponding to liq_alt
    p.create_dataset('cli_g', data=cli_g) # global layer cloud ice water fraction corresponding to alt
    
    p.create_dataset('cl_so', data=cl_so) # southern ocean layer total cloud fraction corresponding to alt
    p.create_dataset('clw_so', data=clw_so) # southern ocean layer cloud liquid water fraction corresponding to liq_alt
    p.create_dataset('cli_so', data=cli_so) # southern ocean layer cloud ice water fraction corresponding to alt
    
    p.create_dataset('lw_frac_alt_lat', data= np.transpose( lw_frac_alt_lat )) # temperature corresponding to alt and lat
    
    p.create_dataset('ta_alt_lat', data= np.transpose( ta_alt_lat )) # temperature corresponding to alt and lat
    p.create_dataset('ta_liq_alt_lat', data= np.transpose( ta_liq_alt_lat )) # temperature corresponding to liq_alt and lat
    p.create_dataset('cl_alt_lat', data=np.transpose( cl_alt_lat ) ) # total cloud fraction corresponding to alt and lat
    p.create_dataset('clw_alt_lat', data=np.transpose( clw_alt_lat ) ) # cloud liquid water fraction corresponding to liq_alt and lat
    p.create_dataset('cli_alt_lat', data=np.transpose( cli_alt_lat ) ) # cloud ice water fraction corresponding to alt and lat

    p.close()

end = time.time()
print('Create dataset took:', end - start, 's')



































