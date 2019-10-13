# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon

This will extract the cloud fraction profile against altitude. 
The global cloud fraction is already averaged over all longitude and latitude at each altitude layer.
[:,0] = altitude
[:,1] = cloud fraction
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
from scipy import ndimage as nd

###############################################################################
start = time.time()

#set location

location = constants.home

class DataSet:
    def __init__( self, type_name ):
        self.type_name = type_name
        self.new_zon_data = np.zeros( constants.lat.size )
        self.new_reg_data = np.zeros( ( constants.lat.size, constants.lon.size ) )
#        self.new_radiation_data = np.zeros( ( constants.lat.size, constants.lon.size ) )
        self.zon_data_counts = np.zeros( constants.lat.size )
        self.reg_data_counts = np.zeros( ( constants.lat.size, constants.lon.size ) )
#        self.radiation_data_counts = np.zeros( ( constants.lat.size, constants.lon.size ) )
        
zon_data_sets = [
        DataSet( 'cld_amount_zon' ),    # divide by 100, alt[4], lat
        DataSet( 'cld_lwp_zon' ),
        DataSet( 'cld_iwp_zon' )
]        

reg_data_sets = [
        DataSet( 'cld_amount_reg' ),    # divide by 100, alt[4], lat
        DataSet( 'cld_lwp_reg' ),
        DataSet( 'cld_iwp_reg' ),  
]      

       
#radiation_data_sets = [
#        DataSet( 'clr_toa_sw_reg' ),       # lat and lon  
#        DataSet( 'clr_toa_lw_reg' ),       # lat and lon  
#        DataSet( 'clr_toa_net_reg' ),       # lat and lon  
#        DataSet( 'all_toa_sw_reg' ),       # lat and lon  
#        DataSet( 'all_toa_lw_reg' ),       # lat and lon  
#        DataSet( 'all_toa_net_reg' )       # lat and lon  
#]        

# The directory where your HDF files are stored
os.chdir(location + 'Data/CERES/Run')

def create_southern_ocean_data( lat, global_data ):
    
    start_index = 0
    end_index = 0
    for index, l in enumerate( lat ):
        if l < -70:
            start_index = index + 1
        if l <= -50:
            end_index = index
    so = global_data[start_index:end_index]
    return so

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
    pprint ( filename )
    f = SD.SD(filename)
    raw_lon = f.select( 'longitude' ).get() + 180
    raw_lat = f.select( 'latitude' ).get()
    raw_lat = np.flip( raw_lat, axis = 0 )
    
    
    for data_set in zon_data_sets:  
        sel = f.select( data_set.type_name )
        data = np.hstack( sel.get()[-1:] )
        fill_value = sel.attributes()['_FillValue']
        
        for l_index, l in enumerate( raw_lat ):
            if l <= constants.min_lat or l >= constants.max_lat:
                continue
            lat_bin = int( ( l - constants.min_lat ) / constants.lat_division)
            val = data[ l_index ]
            if val == fill_value:
                continue
            data_set.new_zon_data[ lat_bin ] += val
            data_set.zon_data_counts[ lat_bin ] += 1


    for data_set in reg_data_sets:  
        sel = f.select( data_set.type_name )
        data = np.hstack(( sel.get()[-1:] ))
        fill_value = sel.attributes()['_FillValue']
        
        for l_index, l in enumerate( raw_lat ):
            if l <= constants.min_lat or l >= constants.max_lat:
                continue
            lat_bin = int( ( l - constants.min_lat ) / constants.lat_division)
            for a_index, a in enumerate( raw_lon ):
                if a < constants.min_lon or a > constants.max_lon:
                    continue
                lon_bin = int( ( a - constants.min_lon ) / constants.lon_division )
                #print( lat_bin, lon_bin )
                val = data[l_index, a_index]
                if val == fill_value:
                    continue
                data_set.new_reg_data[ lat_bin, lon_bin ] += val
                data_set.reg_data_counts[ lat_bin, lon_bin ] += 1
                
                
#    for data_set in radiation_data_sets:  
#        sel = f.select( data_set.type_name )
#        data = sel.get()
#        fill_value = sel.attributes()['_FillValue']
#        
#        for l_index, l in enumerate( raw_lat ):
#            if l <= constants.min_lat or l >= constants.max_lat:
#                continue
#            lat_bin = int( ( l - constants.min_lat ) / constants.lat_division)
#            for a_index, a in enumerate( raw_lon ):
#                if a < constants.min_lon or a > constants.max_lon:
#                    continue
#                lon_bin = int( ( a - constants.min_lon ) / constants.lon_division )
#                #print( lat_bin, lon_bin )
#                val = data[l_index, a_index]
#                if val == fill_value:
#                    continue
#                data_set.new_radiation_data[ lat_bin, lon_bin ] += val
#                data_set.radiation_data_counts[ lat_bin, lon_bin ] += 1

            
for data_set in zon_data_sets:
    data_set.new_zon_data /= data_set.zon_data_counts

for data_set in reg_data_sets:
    data_set.new_reg_data /= data_set.reg_data_counts

#for data_set in radiation_data_sets:
#    data_set.new_radiation_data /= data_set.radiation_data_counts
  
clt = np.flip( fill( zon_data_sets[0].new_zon_data / 100 ), axis = 0 )
clwvi = np.flip( fill( zon_data_sets[1].new_zon_data ), axis = 0 )
clivi = np.flip( fill( zon_data_sets[2].new_zon_data ), axis = 0 )

clt_lat_lon = np.flip( fill( reg_data_sets[0].new_reg_data / 100 ), axis = 0 )
clwvi_lat_lon = np.flip( fill( reg_data_sets[1].new_reg_data ), axis = 0 )
clivi_lat_lon = np.flip( fill( reg_data_sets[2].new_reg_data ), axis = 0 )

#clr_toa_sw_reg = np.flip( fill( radiation_data_sets[0].new_reg_data ), axis = 0 )
#clr_toa_lw_reg = np.flip( fill( radiation_data_sets[1].new_reg_data ), axis = 0 )
#clr_toa_net_reg = np.flip( fill( radiation_data_sets[2].new_reg_data ), axis = 0 )
#all_toa_sw_reg = np.flip( fill( radiation_data_sets[3].new_reg_data ), axis = 0 )
#all_toa_lw_reg = np.flip( fill( radiation_data_sets[4].new_reg_data ), axis = 0 )
#all_toa_net_reg = np.flip( fill( radiation_data_sets[5].new_reg_data ), axis = 0 )

#---create liquid and ice fractions---#

lw_frac = (clwvi / (clwvi + clivi)) * clt
clivi = (clivi / (clwvi + clivi)) * clt
clwvi = lw_frac

lw_lat_lon_frac = (clwvi_lat_lon / (clwvi_lat_lon + clivi_lat_lon)) * clt_lat_lon
clivi_lat_lon = (clivi_lat_lon / (clwvi_lat_lon + clivi_lat_lon)) * clt_lat_lon
clwvi_lat_lon = lw_lat_lon_frac

#
#clr_toa_sw_reg_so = create_southern_ocean_data( constants.lat, clr_toa_sw_reg)
#clr_toa_lw_reg_so = create_southern_ocean_data( constants.lat, clr_toa_lw_reg)
#clr_toa_net_reg_so = create_southern_ocean_data( constants.lat, clr_toa_net_reg)
#all_toa_sw_reg_so = create_southern_ocean_data( constants.lat, all_toa_sw_reg)
#all_toa_lw_reg_so = create_southern_ocean_data( constants.lat, all_toa_lw_reg)
#all_toa_net_reg_so = create_southern_ocean_data( constants.lat, all_toa_net_reg)

#----------------------------#

os.chdir( location + '/climate-analysis/reduced_data' )

save_filename = 'Jun_2006_Jun_2011_CERES.h5'

with h5py.File(save_filename, 'w') as p:
            
    p.create_dataset('clt', data=clt)
    p.create_dataset('clwvi', data=clwvi) 
    p.create_dataset('clivi', data=clivi) 
    
    p.create_dataset('clt_lat_lon', data=clt_lat_lon)
    p.create_dataset('clwvi_lat_lon', data=clwvi_lat_lon) 
    p.create_dataset('clivi_lat_lon', data=clivi_lat_lon) 

#    p.create_dataset('clr_toa_sw_reg', data=clr_toa_sw_reg) 
#    p.create_dataset('clr_toa_sw_reg_so', data=clr_toa_sw_reg_so) 
#    p.create_dataset('clr_toa_lw_reg', data=clr_toa_lw_reg) 
#    p.create_dataset('clr_toa_lw_reg_so', data=clr_toa_lw_reg_so) 
#    p.create_dataset('clr_toa_net_reg', data=clr_toa_net_reg) 
#    p.create_dataset('clr_toa_net_reg_so', data=clr_toa_net_reg_so) 
#    p.create_dataset('all_toa_sw_reg', data=all_toa_sw_reg) 
#    p.create_dataset('all_toa_sw_reg_so', data=all_toa_sw_reg_so) 
#    p.create_dataset('all_toa_lw_reg', data=all_toa_lw_reg) 
#    p.create_dataset('all_toa_lw_reg_so', data=all_toa_lw_reg_so) 
#    p.create_dataset('all_toa_net_reg', data=all_toa_net_reg) 
#    p.create_dataset('all_toa_net_reg_so', data=all_toa_net_reg_so) 

    p.close()

end = time.time()
print('Create dataset took:', end - start, 's')

