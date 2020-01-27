# -*- coding: utf-8 -*-
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

#read hdf5
import h5py
import os
import matplotlib.pyplot as plt
import constants
import numpy as np
import cartopy.crs as ccrs

#--- Set Location, date period and model ---#

# specify model from the list above
# specify location: home, uni, hdd, laptop
location = constants.uni + 'climate-analysis/reduced_data'

os.chdir( location )


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


h5f = h5py.File( constants.date_ceres + '_CERES' + '.h5', 'r')
cer_albedo_reg = h5f['albedo_reg'][:]
cer_rtmt_lat_lon = h5f['rtmt_lat_lon'][:]
cer_all_toa_net_so = h5f['all_toa_net_so'][:]
cer_albedo_so = h5f['albedo_so'][:]

cer_net_glob = np.nanmean( np.nanmean( cer_rtmt_lat_lon, axis = 0 ), axis = 0 )
cer_albedo_glob = np.nanmean( np.nanmean( cer_albedo_reg, axis = 0 ), axis = 0 )
cer_albedo_so = np.nanmean( cer_albedo_so, axis = 0 )
cer_all_toa_net_so = np.nanmean( cer_all_toa_net_so, axis = 0 )

print( 'CERES Net global energy = ' + str(cer_net_glob) )
print( 'CERES Net SO energy = ' + str(cer_all_toa_net_so) )
print( 'CERES Global albedo = ' + str(cer_albedo_glob) )
print( 'CERES SO albedo = ' + str(cer_albedo_so) )


for model in constants.all_models:

    if model == 'CMIP5-CESM1-CAM5' or model == 'CMIP5-GFDL-HIRAM-C360' or model == 'CMIP5-GISS-E2R' or model == 'CMIP5-IPSL-CM5A-LR' or model == 'CMIP5-MIROC5' or model == 'CMIP5-MRI-CGCM3': 
        h5f = h5py.File( constants.date_cmip5 + '_' + model + '.h5', 'r')
        rtmt_lat_lon = h5f['rtmt_lat_lon'][:]
        albedo_reg = h5f['albedo_reg'][:]
        albedo_so = h5f['albedo_so'][:]

    else:
        h5f = h5py.File( constants.date_cmip6 + '_' + model + '.h5', 'r')
        rtmt_lat_lon = h5f['rtmt_lat_lon'][:]
        albedo_reg = h5f['albedo_reg'][:]
        albedo_so = h5f['albedo_so'][:]

    toa_net_so = np.nanmean( np.nanmean( create_southern_ocean_data( constants.lat, rtmt_lat_lon ), axis = 0 ), axis = 0 )


    net_glob = np.nanmean( np.nanmean( rtmt_lat_lon, axis = 0 ), axis = 0 )
    albedo_glob = np.nanmean( np.nanmean( albedo_reg, axis = 0 ), axis = 0 )
    albedo_so = np.nanmean( albedo_so, axis = 0 )

    net_glob_bias = cer_net_glob - net_glob
    net_so_bias = cer_all_toa_net_so - toa_net_so
    albedo_glob_bias = cer_albedo_glob - albedo_glob
    albedo_so_bias = cer_albedo_so - albedo_so

    print( model )
    print( 'Net global energy = ' + str(net_glob) )
    print( 'Global albedo = ' + str(albedo_glob) )
    print( 'SO albedo = ' + str(albedo_so) )


    print( 'Global energy bias = ' + str(net_glob_bias))
    print( 'SO energy bias = ' + str(net_so_bias) )
    print( 'global albedo bias = ' + str(albedo_glob_bias) )
    print( 'SO albedo bias = ' + str(albedo_so_bias) )



# ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
# ax.coastlines()
# p = ax.contourf(constants.lon, constants.lat, cer_albedo_reg, transform=ccrs.PlateCarree(), cmap='coolwarm')
# cbar = plt.colorbar(p, orientation='horizontal')
# cbar.set_label('Cloud Fraction')
# ax.set_title('Total Cloud Fraction - CMIP6-GFDL-AM4')

# plt.show()




