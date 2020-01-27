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

lc_zon_data_sets = [
        DataSet( 'cld_amount_zon' ),    # divide by 100, alt[3], lat
        DataSet( 'cld_lwp_zon' ),
        DataSet( 'cld_iwp_zon' )
]        

lc_reg_data_sets = [
        DataSet( 'cld_amount_reg' ),    # divide by 100, alt[3], lat
        DataSet( 'cld_lwp_reg' ),
        DataSet( 'cld_iwp_reg' ),  
]      

       
reg_radiation_data_sets = [
        DataSet( 'clr_toa_sw_reg' ),       # lat and lon  
        DataSet( 'clr_toa_lw_reg' ),       # lat and lon  
        DataSet( 'clr_toa_net_reg' ),       # lat and lon  
        DataSet( 'all_toa_sw_reg' ),       # lat and lon  
        DataSet( 'all_toa_lw_reg' ),       # lat and lon  
        DataSet( 'all_toa_net_reg' ),       # lat and lon  
        DataSet( 'toa_sw_insol_reg' ),       # lat and lon  
        DataSet( 'all_toa_alb_reg' ),       # lat and lon 
]        

zon_radiation_data_sets = [
        DataSet( 'clr_toa_sw_zon' ),       # lat 
        DataSet( 'clr_toa_lw_zon' ),       # lat
        DataSet( 'clr_toa_net_zon' ),       # lat
        DataSet( 'all_toa_sw_zon' ),       # lat
        DataSet( 'all_toa_lw_zon' ),       # lat
        DataSet( 'all_toa_net_zon' ),       # lat
        DataSet( 'toa_sw_insol_zon' ),       # lat
        DataSet( 'all_toa_alb_zon' ),       # lat
]        

glob_radiation_data_sets = [
        DataSet( 'clr_toa_sw_glob' ), 
        DataSet( 'clr_toa_lw_glob' ), 
        DataSet( 'clr_toa_net_glob' ),
        DataSet( 'all_toa_sw_glob' ),
        DataSet( 'all_toa_lw_glob' ),
        DataSet( 'all_toa_net_glob' ), 
        DataSet( 'toa_sw_insol_glob' ),
        DataSet( 'all_toa_alb_glob' ),
]        

# The directory where your HDF files are stored
os.chdir(location + 'Data/CERES/Run' )

# Load every file in the directory

#Load the first file, and use it to get the longditude and latitude data size
files = os.listdir()
reference_file = SD.SD( files[ 0 ] )
raw_lon = reference_file.select( 'longitude' ).get()
raw_lat = reference_file.select( 'latitude' ).get()

for data_set in zon_data_sets:
    data_set.new_zon_data = np.zeros( raw_lat.size )
    data_set.zon_data_counts = np.zeros( raw_lat.size )
for data_set in lc_zon_data_sets:
    data_set.new_lc_zon_data = np.zeros( raw_lat.size )
    data_set.lc_zon_data_counts = np.zeros( raw_lat.size )
for data_set in reg_data_sets:
    data_set.new_reg_data = np.zeros( ( raw_lat.size, raw_lon.size ) )
    data_set.reg_data_counts = np.zeros( ( raw_lat.size, raw_lon.size ) )
for data_set in lc_reg_data_sets:
    data_set.new_lc_reg_data = np.zeros( ( raw_lat.size, raw_lon.size ) )
    data_set.lc_reg_data_counts = np.zeros( ( raw_lat.size, raw_lon.size ) )
for data_set in reg_radiation_data_sets:
    data_set.new_reg_radiation_data = np.zeros( ( raw_lat.size, raw_lon.size ) )
    data_set.reg_radiation_data_counts = np.zeros( ( raw_lat.size, raw_lon.size ) )
for data_set in zon_radiation_data_sets:
    data_set.new_zon_radiation_data = np.zeros( ( raw_lat.size ) )
    data_set.zon_radiation_data_counts = np.zeros( ( raw_lat.size ) )
for data_set in glob_radiation_data_sets:
    data_set.new_glob_radiation_data = 0
    data_set.glob_radiation_data_counts = 0


for filename in os.listdir():
    if os.path.isdir( filename ):
        continue
    pprint ( filename )
    f = SD.SD(filename)
    raw_lon = f.select( 'longitude' ).get() + 180
    raw_lat = f.select( 'latitude' ).get()
    raw_lat = np.flip( raw_lat, axis = 0 )
    idx = np.where(raw_lon == 180)
    start_idx = np.abs(raw_lat - (-69.5)).argmin()
    end_idx = np.abs(raw_lat - (-50)).argmin()


    for data_set in zon_data_sets:  

        sel = f.select( data_set.type_name )
        data = sel.get()[4]

        fill_value = sel.attributes()['_FillValue']
        
        for l_index, l in enumerate( raw_lat ):
            val = data[ l_index ]
            if val == fill_value:
                continue
            data_set.new_zon_data[ l_index ] += val
            data_set.zon_data_counts[ l_index ] += 1


    for data_set in lc_zon_data_sets:  

        sel = f.select( data_set.type_name )
        data = sel.get()[3]
        fill_value = sel.attributes()['_FillValue']
        
        for l_index, l in enumerate( raw_lat ):
            val = data[ l_index ]
            if val == fill_value:
                continue
            data_set.new_lc_zon_data[ l_index ] += val
            data_set.lc_zon_data_counts[ l_index ] += 1


    for data_set in reg_data_sets:  

        sel = f.select( data_set.type_name )
        data = sel.get()[4]
        data[:] = np.roll(data[:], 179)
        fill_value = sel.attributes()['_FillValue']
        
        for l_index, l in enumerate( raw_lat ):
            for a_index, a in enumerate( raw_lon ):
                val = data[l_index, a_index]
                if val == fill_value:
                    continue
                data_set.new_reg_data[ l_index, a_index ] += val
                data_set.reg_data_counts[ l_index, a_index ] += 1
                

    for data_set in lc_reg_data_sets:  

        sel = f.select( data_set.type_name )
        data = sel.get()[3]
        data[:] = np.roll(data[:], 179)
        fill_value = sel.attributes()['_FillValue']
        
        for l_index, l in enumerate( raw_lat ):
            for a_index, a in enumerate( raw_lon ):
                val = data[l_index, a_index]
                if val == fill_value:
                    continue
                data_set.new_lc_reg_data[ l_index, a_index ] += val
                data_set.lc_reg_data_counts[ l_index, a_index ] += 1


    for data_set in reg_radiation_data_sets:
        sel = f.select( data_set.type_name )
        data = sel.get()[:]
        data[:] = np.roll(data[:], 179)
        fill_value = sel.attributes()['_FillValue']

        for l_index, l in enumerate( raw_lat ):
            for a_index, a in enumerate( raw_lon ):
                val = data[l_index, a_index]
                if val == fill_value:
                    continue
                data_set.new_reg_radiation_data[ l_index, a_index ] += val
                data_set.reg_radiation_data_counts[ l_index, a_index ] += 1


    for data_set in zon_radiation_data_sets:
        sel = f.select( data_set.type_name )
        data = sel.get()[:]
        fill_value = sel.attributes()['_FillValue']

        for l_index, l in enumerate( raw_lat ):
                val = data[l_index]
                if val == fill_value:
                    continue
                data_set.new_zon_radiation_data[ l_index ] += val
                data_set.zon_radiation_data_counts[ l_index ] += 1


    for data_set in glob_radiation_data_sets:
        sel = f.select( data_set.type_name )
        data = sel.get()[:]
        fill_value = sel.attributes()['_FillValue']
        if data == fill_value:
            continue
        data_set.new_glob_radiation_data += data
        data_set.glob_radiation_data_counts += 1

            
for data_set in zon_data_sets:
    data_set.new_zon_data /= data_set.zon_data_counts

for data_set in reg_data_sets:
    data_set.new_reg_data /= data_set.reg_data_counts

for data_set in lc_zon_data_sets:
    data_set.new_lc_zon_data /= data_set.lc_zon_data_counts

for data_set in lc_reg_data_sets:
    data_set.new_lc_reg_data /= data_set.lc_reg_data_counts

for data_set in reg_radiation_data_sets:
    data_set.new_reg_radiation_data /= data_set.reg_radiation_data_counts

for data_set in zon_radiation_data_sets:
    data_set.new_zon_radiation_data /= data_set.zon_radiation_data_counts

for data_set in glob_radiation_data_sets:
    data_set.new_glob_radiation_data /= data_set.glob_radiation_data_counts

clt = np.flip( zon_data_sets[0].new_zon_data, axis = 0 ) / 100
clwvi = np.flip( zon_data_sets[1].new_zon_data, axis = 0 ) / 1000
clivi = np.flip( zon_data_sets[2].new_zon_data, axis = 0 ) / 1000

clt_lc = np.flip( lc_zon_data_sets[0].new_lc_zon_data, axis = 0 ) / 100
clwvi_lc = np.flip( lc_zon_data_sets[1].new_lc_zon_data, axis = 0 ) / 1000
clivi_lc = np.flip( lc_zon_data_sets[2].new_lc_zon_data, axis = 0 ) / 1000

clt_lc_lat_lon = np.flip( lc_reg_data_sets[0].new_lc_reg_data, axis = 0 ) / 100
clwvi_lc_lat_lon = np.flip( lc_reg_data_sets[1].new_lc_reg_data, axis = 0 ) / 1000
clivi_lc_lat_lon = np.flip( lc_reg_data_sets[2].new_lc_reg_data, axis = 0 ) / 1000

clt_lat_lon = np.flip( reg_data_sets[0].new_reg_data, axis = 0 ) / 100
clwvi_lat_lon = np.flip( reg_data_sets[1].new_reg_data, axis = 0 ) / 1000
clivi_lat_lon = np.flip( reg_data_sets[2].new_reg_data, axis = 0 ) / 1000

clr_toa_sw_reg = np.flip( reg_radiation_data_sets[0].new_reg_radiation_data, axis = 0 )
clr_toa_lw_reg = np.flip( reg_radiation_data_sets[1].new_reg_radiation_data, axis = 0 )
clr_toa_net_reg = np.flip( reg_radiation_data_sets[2].new_reg_radiation_data, axis = 0 )
all_toa_sw_reg = np.flip( reg_radiation_data_sets[3].new_reg_radiation_data, axis = 0 )
all_toa_lw_reg = np.flip( reg_radiation_data_sets[4].new_reg_radiation_data, axis = 0 )
all_toa_net_reg = np.flip( reg_radiation_data_sets[5].new_reg_radiation_data, axis = 0 )
all_toa_sw_insol_reg = np.flip( reg_radiation_data_sets[6].new_reg_radiation_data, axis = 0 )
all_toa_alb_reg = np.flip( reg_radiation_data_sets[7].new_reg_radiation_data, axis = 0 )

clr_toa_sw_zon= np.flip( zon_radiation_data_sets[0].new_zon_radiation_data, axis = 0 )
clr_toa_lw_zon= np.flip( zon_radiation_data_sets[1].new_zon_radiation_data, axis = 0 )
clr_toa_net_zon= np.flip( zon_radiation_data_sets[2].new_zon_radiation_data, axis = 0 )
all_toa_sw_zon= np.flip( zon_radiation_data_sets[3].new_zon_radiation_data, axis = 0 )
all_toa_lw_zon= np.flip( zon_radiation_data_sets[4].new_zon_radiation_data, axis = 0 )
all_toa_net_zon= np.flip( zon_radiation_data_sets[5].new_zon_radiation_data, axis = 0 )
all_toa_sw_insol_zon= np.flip( zon_radiation_data_sets[6].new_zon_radiation_data, axis = 0 )
all_toa_alb_zon= np.flip( zon_radiation_data_sets[7].new_zon_radiation_data, axis = 0 )

clr_toa_sw_glob = glob_radiation_data_sets[0].new_glob_radiation_data
clr_toa_lw_glob = glob_radiation_data_sets[1].new_glob_radiation_data
clr_toa_net_glob = glob_radiation_data_sets[2].new_glob_radiation_data
all_toa_sw_glob = glob_radiation_data_sets[3].new_glob_radiation_data
all_toa_lw_glob = glob_radiation_data_sets[4].new_glob_radiation_data
all_toa_net_glob = glob_radiation_data_sets[5].new_glob_radiation_data
all_toa_sw_insol_glob = glob_radiation_data_sets[6].new_glob_radiation_data
albedo_glob = glob_radiation_data_sets[7].new_glob_radiation_data

#------#

clr_toa_sw_so = np.mean( constants.create_southern_ocean_data( raw_lat, clr_toa_sw_zon), axis = 0 )
clr_toa_lw_so = np.mean( constants.create_southern_ocean_data( raw_lat, clr_toa_lw_zon), axis = 0 )
clr_toa_net_so = np.mean( constants.create_southern_ocean_data( raw_lat, clr_toa_net_zon), axis = 0 )
all_toa_sw_so_zon = constants.create_southern_ocean_data( raw_lat, all_toa_sw_zon)
all_toa_lw_so_zon = constants.create_southern_ocean_data( raw_lat, all_toa_lw_zon)
all_toa_net_so_zon = constants.create_southern_ocean_data( raw_lat, all_toa_net_zon)
all_toa_sw_insol_so_zon = constants.create_southern_ocean_data( raw_lat, all_toa_sw_insol_zon)
albedo_so_zon = constants.create_southern_ocean_data( raw_lat, all_toa_alb_zon)

cre_sw_reg = clr_toa_sw_reg - all_toa_sw_reg
cre_lw_reg = clr_toa_lw_reg - all_toa_lw_reg

cre_sw = clr_toa_sw_glob - all_toa_sw_glob
cre_lw = clr_toa_lw_glob - all_toa_lw_glob

cre_sw_reg_so = constants.create_southern_ocean_data( raw_lat, cre_sw_reg )
cre_lw_reg_so = constants.create_southern_ocean_data( raw_lat, cre_lw_reg )
albedo_so_reg = constants.create_southern_ocean_data( raw_lat, all_toa_alb_reg )
rtmt_so_reg = constants.create_southern_ocean_data( raw_lat, all_toa_net_reg )


print( 'Incoming Solar Energy = ' + str( all_toa_sw_insol_glob ) )
print( 'Outgoing SW Energy = ' + str( all_toa_sw_glob ) )
print( 'Outgoing LW Energy = ' + str( all_toa_lw_glob ) )
print( 'Net Global Energy Flux = ' + str( all_toa_net_glob ) )
print( 'Net SO Energy Flux = ' + str( constants.global2DMean(rtmt_so_reg, raw_lat[start_idx:end_idx]) ) )
print( 'Global Albedo = ' + str( albedo_glob ) )
print( 'SO albedo = ' + str( constants.global2DMean(albedo_so_reg, raw_lat[start_idx:end_idx]) ) )

print( 'Global CRE SW = ' + str( cre_sw ) )
print( 'Global CRE LW = ' + str( cre_lw ) )
print( 'SO CRE SW = ' + str( constants.global2DMean(cre_sw_reg_so, raw_lat[start_idx:end_idx]) ) )
print( 'SO CRE LW = ' + str( constants.global2DMean(cre_lw_reg_so, raw_lat[start_idx:end_idx]) ) )



#----------------------------#

os.chdir( location + 'climate-analysis/reduced_data' )

save_filename = 'all_CERES.h5'

with h5py.File(save_filename, 'w') as p:
            
    p.create_dataset('clt', data=clt)
    p.create_dataset('clwvi', data=clwvi) 
    p.create_dataset('clivi', data=clivi) 
    
    p.create_dataset('clt_lc', data=clt_lc)
    p.create_dataset('clwvi_lc', data=clwvi_lc) 
    p.create_dataset('clivi_lc', data=clivi_lc) 

    p.create_dataset('clt_lat_lon', data=clt_lat_lon)
    p.create_dataset('clwvi_lat_lon', data=clwvi_lat_lon) 
    p.create_dataset('clivi_lat_lon', data=clivi_lat_lon) 

    p.create_dataset('clt_lc_lat_lon', data=clt_lc_lat_lon)
    p.create_dataset('clwvi_lc_lat_lon', data=clwvi_lc_lat_lon) 
    p.create_dataset('clivi_lc_lat_lon', data=clivi_lc_lat_lon) 

    p.create_dataset('albedo_reg', data=all_toa_alb_reg) 
    p.create_dataset('rtmt_lat_lon', data=all_toa_net_reg) 

    p.create_dataset('all_toa_net_so', data=all_toa_net_so_zon) 

    p.create_dataset('albedo_so', data=albedo_so_zon) 

    p.close()

end = time.time()
print('Create dataset took:', end - start, 's')

