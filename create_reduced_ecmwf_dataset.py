# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon


"""
import numpy as np
import os
from netCDF4 import Dataset
from netCDF4 import date2index
import h5py
import math
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline
import datetime
from pprint import pprint
from sklearn.impute import SimpleImputer
import constants
import time
from scipy import ndimage as nd

###############################################################################
start = time.time()

#set location

location = constants.home
os.chdir( location + 'Data/ECMWF/pressure_levels' )

#lat = constants.lat


class DataSets:
    def __init__( self, type_name ):
        self.type_name = type_name
        self.new_data = np.zeros(( constants.lat.size, constants.alt.size ))
        self.data_counts = np.zeros( (constants.lat.size, constants.alt.size ))
        
data_sets = [
        DataSets('cc'),
        DataSets('clwc'),
        DataSets('ciwc'),
        DataSets('t')  ,
]        

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

def extract_data( type, f ):
    return np.array( f.variables[type][:] )


def fill(data, invalid=None):
    """
    Replace the value of invalid 'data' cells (indicated by 'invalid') 
    by the value of the nearest valid data cell

    Input:
        data:    numpy array of any dimension
        invalid: a binary array of same shape as 'data'. 
                 data value are replaced where invalid is True
                 If None (default), use: invalid  = np.isnan(data)

    Output: 
        Return a filled array. 
    """    
    if invalid is None: invalid = np.isnan(data)

    ind = nd.distance_transform_edt(invalid, 
                                    return_distances=False, 
                                    return_indices=True)
    return data[tuple(ind)]


# Load every file in the directory
for filename in os.listdir(): 
    with Dataset( filename, 'r') as f:
        p = np.flip( extract_data( 'level', f ) )        
        raw_alt = np.empty((p.size,1),dtype=float)
        state = 0
        i = 0
        for item in p:
            if state == 0:
                newalt = (288.19 - 288.08*((item/1012.90)**(1/5.256)))/6.49
                if newalt > 11:
                    state = 1
            if state == 1:
                newalt = (1.73 - math.log(item/226.50))/0.157
                if( newalt > 25 ):
                    state = 2
            if state == 2:
                newalt = (216.6*((item/24.88)**(1/-11.388)) - 141.94)/2.99
            raw_alt[i] = newalt
            i+=1
        raw_alt = np.transpose( raw_alt )[0]
        raw_lat = np.flip( extract_data( 'latitude', f ) )
        for data_set in data_sets:  
            data = extract_data( data_set.type_name, f )
            data = np.mean( data, axis = 0 )
            data = np.transpose( np.flip( np.mean( data, axis = -1 ) ) )
            
            for l_index, l in enumerate( raw_lat ):
                if l <= constants.min_lat or l >= constants.max_lat:
                    continue
                lat_bin = int( ( l - constants.min_lat ) / constants.lat_division)
                for a_index, a in enumerate( raw_alt ):
                    if a < constants.min_alt or a > constants.max_alt:
                        continue
                    alt_bin = int( ( a - constants.min_alt ) / constants.alt_division )
                    #print( lat_bin, alt_bin )
                    val = data[l_index, a_index]
                    data_set.new_data[ lat_bin, alt_bin ] += val
                    data_set.data_counts[ lat_bin, alt_bin ] += 1
            
for data_set in data_sets:
    data_set.new_data /= data_set.data_counts
  
cl_alt_lat = fill( data_sets[0].new_data )
clw_alt_lat = fill( data_sets[1].new_data )
cli_alt_lat = fill( data_sets[2].new_data )
ta_alt_lat = fill( data_sets[3].new_data )

#---create liquid and ice fractions---#

lw_frac_alt_lat = (clw_alt_lat / (clw_alt_lat + cli_alt_lat)) * cl_alt_lat
cli_alt_lat = (cli_alt_lat / (clw_alt_lat + cli_alt_lat)) * cl_alt_lat

#---create reduced altitude liquid fraction and temperatures---#

interpolated = interpolate.interp2d( constants.alt, constants.lat, lw_frac_alt_lat, kind = 'cubic')
clw_alt_lat = interpolated( constants.liq_alt, constants.lat )

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

save_filename = 'ECMWF.h5'

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

