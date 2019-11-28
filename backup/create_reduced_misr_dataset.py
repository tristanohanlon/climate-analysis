# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon - University of Auckland & Jonathan Rogers


"""

import time
import numpy as np
import os
from pyhdf import SD
import h5py
from scipy import integrate
from sklearn.impute import SimpleImputer
import constants
from pprint import pprint


#set location

location = constants.home + '/Data/MISR/MIL3YCFA/'
os.chdir( location )
start_year = 2006 # inclusive
end_year = 2011 # inclusive

raw_lat = np.arange( -89.75, 90, 0.5 )
raw_alt = np.arange( -0.750, 21.750, 0.50 )

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


############################################################################### Get cl (lat, lon, alt)

new_data = np.zeros(( constants.lat.size, constants.alt.size ))
data_counts = np.zeros( (constants.lat.size, constants.alt.size ))


for file in os.listdir():
    f = SD.SD(file)
    sel = f.select('CorrCloudTopHeightFraction_Avg')

    data = sel.get()
    data[data < 0] = np.nan # convert all fill values to nan
    data = np.nanmean( data, axis = 1 ) # average over longitude
    
    for l_index, l in enumerate( raw_lat ):
        if l < constants.min_lat or l > constants.max_lat:
            continue
        lat_bin = int( ( l - constants.min_lat ) / constants.lat_division)
        for a_index, a in enumerate( raw_alt ):
            if a < constants.min_alt or a > constants.max_alt:
                continue
            alt_bin = int( ( a - constants.min_alt ) / constants.alt_division )
            val = data[l_index, a_index]
            if np.isnan(val):
                continue
            new_data[ lat_bin, alt_bin ] += val
            data_counts[ lat_bin, alt_bin ] += 1
            
new_data /= data_counts

############################################################################### create global and southern ocean variables

cl_so = create_southern_ocean_data( constants.lat, new_data)
cl_g = np.nanmean( new_data, axis = 0 )

os.chdir( constants.home + '/climate-analysis/reduced_data' )

save_filename = 'MISR.h5'

with h5py.File(save_filename, 'w') as p:
    p.create_dataset('cl_g', data=cl_g) # global layer total cloud fraction corresponding to alt
    p.create_dataset('cl_so', data=cl_so) # southern ocean layer total cloud fraction corresponding to alt
    p.close()


