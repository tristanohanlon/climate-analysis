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
model = 'CMIP6-MIROC6' 

# specify location: home, uni, hdd, laptop
location = constants.home + 'climate-analysis/reduced_data'

os.chdir( location )

if model == 'CCCM':
    h5f = h5py.File( '1' + constants.date_cccm + '_' + model + '.h5', 'r')
    cl_alt_lat = h5f['cl_alt_lat'][:]
    cl_g = h5f['cl_g'][:]
    cl_so = h5f['cl_so'][:]
    cli_alt_lat = h5f['cli_alt_lat'][:]
    cli_g = h5f['cli_g'][:]
    cli_so = h5f['cli_so'][:]
    clw_alt_lat = h5f['clw_alt_lat'][:]
    clw_g = h5f['clw_g'][:]
    clw_so = h5f['clw_so'][:]
    ta_alt_lat = h5f['ta_alt_lat'][:]
    clw_t_g = h5f['clw_t_g'][:]
    clw_t_so = h5f['clw_t_so'][:]
    full_clw_alt_lat = h5f['full_clw_alt_lat'][:]
    full_ta_alt_lat = h5f['full_ta_alt_lat'][:]
    clt = h5f['clt'][:]
    clt_lat_lon = h5f['clt_lat_lon'][:]
    clwvi = h5f['clwvi'][:]
    clwvi_lat_lon = h5f['clwvi_lat_lon'][:]
    clivi = h5f['clivi'][:]
    clivi_lat_lon = h5f['clivi_lat_lon'][:]
   
elif model == 'CERES':
    h5f = h5py.File( constants.date_ceres + '_' + model + '.h5', 'r')
    clt = h5f['clt'][:]
    clt_lat_lon = h5f['clt_lat_lon'][:]
    clwvi = h5f['clwvi'][:]
    clwvi_lat_lon = h5f['clwvi_lat_lon'][:]
    clivi = h5f['clivi'][:]
    clivi_lat_lon = h5f['clivi_lat_lon'][:]

    clt_lc = h5f['clt_lc'][:]
    clt_lc_lat_lon = h5f['clt_lc_lat_lon'][:]
    clwvi_lc = h5f['clwvi_lc'][:]
    clwvi_lc_lat_lon = h5f['clwvi_lc_lat_lon'][:]
    clivi_lc = h5f['clivi_lc'][:]
    clivi_lc_lat_lon = h5f['clivi_lc_lat_lon'][:]

    albedo_reg = h5f['albedo_reg'][:]
    rtmt_lat_lon = h5f['rtmt_lat_lon'][:]

    all_toa_net_so = h5f['all_toa_net_so']

    albedo_so = h5f['albedo_so']

    
elif model == 'ECMWF':    
    h5f = h5py.File( constants.date_cmip6 + '_' + model + '.h5', 'r')
    cl_alt_lat = h5f['cl_alt_lat'][:]
    cl_g = h5f['cl_g'][:]
    cl_so = h5f['cl_so'][:]
    cli_alt_lat = h5f['cli_alt_lat'][:]
    cli_g = h5f['cli_g'][:]
    cli_so = h5f['cli_so'][:]
    clw_alt_lat = h5f['clw_alt_lat'][:]
    clw_g = h5f['clw_g'][:]
    clw_so = h5f['clw_so'][:]
    ta_alt_lat = h5f['ta_alt_lat'][:]
    clw_t_g = h5f['clw_t_g'][:]
    clw_t_so = h5f['clw_t_so'][:]
    full_clw_alt_lat = h5f['full_clw_alt_lat'][:]
    full_ta_alt_lat = h5f['full_ta_alt_lat'][:]
    clt = h5f['clt'][:]
    clt_lat_lon = h5f['clt_lat_lon'][:]
    clwvi = h5f['clwvi'][:]
    clwvi_lat_lon = h5f['clwvi_lat_lon'][:]
    clivi = h5f['clivi'][:]
    clivi_lat_lon = h5f['clivi_lat_lon'][:]
    clt_lc = h5f['clt_lc'][:]
    clt_lc_lat_lon = h5f['clt_lc_lat_lon'][:]
 
    w_g = h5f['w_g'][:]
    w_so = h5f['w_so'][:]
    w_alt_lat = h5f['w_alt_lat'][:]

elif model == 'CALIPSO':    
    h5f = h5py.File( 'Jun_2006_Jun_2011_CALIPSO.h5', 'r')
    cl_alt_lat = h5f['cl_alt_lat'][:]
    cl_g = h5f['cl_g'][:]
    cl_so = h5f['cl_so'][:]
    cli_alt_lat = h5f['cli_alt_lat'][:]
    cli_g = h5f['cli_g'][:]
    cli_so = h5f['cli_so'][:]
    clw_g = h5f['clw_g'][:]
    clw_so = h5f['clw_so'][:]
    clw_t_g = h5f['clw_t_g'][:]
    clw_t_so = h5f['clw_t_so'][:]
    clw_alt_lat = h5f['clw_alt_lat'][:]
    clt = h5f['clt'][:]
    clt_lat_lon = h5f['clt_lat_lon'][:]
    clwvi = h5f['clwvi'][:]
    clwvi_lat_lon = h5f['clwvi_lat_lon'][:]
    clivi = h5f['clivi'][:]
    clivi_lat_lon = h5f['clivi_lat_lon'][:]
    full_clw_alt_lat = h5f['full_clw_alt_lat'][:]
    full_ta_alt_lat = h5f['full_ta_alt_lat'][:]
    ta_alt_lat = h5f['ta_alt_lat'][:]
    clt_lc = h5f['clt_lc'][:]
    clwvi_lc = h5f['clwvi_lc'][:]
    clt_lc_lat_lon = h5f['clt_lc_lat_lon'][:]
    clwvi_lc_lat_lon = h5f['clwvi_lc_lat_lon'][:]

    
elif model == 'CMIP5-CESM1-CAM5' or model == 'CMIP5-GFDL-HIRAM-C360' or model == 'CMIP5-GISS-E2R' or model == 'CMIP5-IPSL-CM5A-LR' or model == 'CMIP5-MIROC5' or model == 'CMIP5-MRI-CGCM3': 
    h5f = h5py.File( constants.date_cmip5 + '_' + model + '.h5', 'r')
    cl_alt_lat = h5f['cl_alt_lat'][:]
    cl_g = h5f['cl_g'][:]
    cl_so = h5f['cl_so'][:]
    cli_alt_lat = h5f['cli_alt_lat'][:]
    cli_g = h5f['cli_g'][:]
    cli_so = h5f['cli_so'][:]
    clt = h5f['clt'][:]
    clw_alt_lat = h5f['clw_alt_lat'][:]
    clw_g = h5f['clw_g'][:]
    clw_so = h5f['clw_so'][:]
    ta_alt_lat = h5f['ta_alt_lat'][:]
    clw_t_g = h5f['clw_t_g'][:]
    clw_t_so = h5f['clw_t_so'][:]
    full_clw_alt_lat = h5f['full_clw_alt_lat'][:]
    full_ta_alt_lat = h5f['full_ta_alt_lat'][:]
    clt = h5f['clt'][:]
    clt_lat_lon = h5f['clt_lat_lon'][:]
    clwvi = h5f['clwvi'][:]
    clwvi_lat_lon = h5f['clwvi_lat_lon'][:]
    clivi = h5f['clivi'][:]
    clivi_lat_lon = h5f['clivi_lat_lon'][:]

    clt_lc = h5f['clt_lc'][:]
    clwvi_lc = h5f['clwvi_lc'][:]
    clt_lc_lat_lon = h5f['clt_lc_lat_lon'][:]
    clwvi_lc_lat_lon = h5f['clwvi_lc_lat_lon'][:]

    rsdt_lat_lon = h5f['rsdt_lat_lon'][:]
    rsut_lat_lon = h5f['rsut_lat_lon'][:]
    rsutcs_lat_lon = h5f['rsutcs_lat_lon'][:]
    rtmt_lat_lon = h5f['rtmt_lat_lon'][:]
    albedo_reg = h5f['albedo_reg'][:]
    albedo_so = h5f['albedo_so']

    if model == 'CMIP5-MRI-CGCM3':
        rlut_lat_lon = h5f['rlut_lat_lon'][:]
        rlutcs_lat_lon = h5f['rlutcs_lat_lon'][:]

    if model == 'CMIP5-GISS-E2R' or model == 'CMIP5-MIROC5':
        mmrdust_lat_lon = h5f['mmrdust_lat_lon'][:]
        mmroa_lat_lon = h5f['mmroa_lat_lon'][:]
        mmrso4_lat_lon = h5f['mmrso4_lat_lon'][:]
        mmrbc_lat_lon = h5f['mmrbc_lat_lon'][:]
        mmrss_lat_lon = h5f['mmrss_lat_lon'][:]
        aerosol_norm_lat_lon = h5f['aerosol_norm_lat_lon'][:]


else:
    h5f = h5py.File( constants.date_cmip6 + '_' + model + '.h5', 'r')
    cl_alt_lat = h5f['cl_alt_lat'][:]
    cl_g = h5f['cl_g'][:]
    cl_so = h5f['cl_so'][:]
    cli_alt_lat = h5f['cli_alt_lat'][:]
    cli_g = h5f['cli_g'][:]
    cli_so = h5f['cli_so'][:]
    clt = h5f['clt'][:]
    clw_alt_lat = h5f['clw_alt_lat'][:]
    clw_g = h5f['clw_g'][:]
    clw_so = h5f['clw_so'][:]
    ta_alt_lat = h5f['ta_alt_lat'][:]
    clw_t_g = h5f['clw_t_g'][:]
    clw_t_so = h5f['clw_t_so'][:]
    full_clw_alt_lat = h5f['full_clw_alt_lat'][:]
    full_ta_alt_lat = h5f['full_ta_alt_lat'][:]
    clt = h5f['clt'][:]
    clt_lat_lon = h5f['clt_lat_lon'][:]
    clwvi = h5f['clwvi'][:]
    clwvi_lat_lon = h5f['clwvi_lat_lon'][:]
    clivi = h5f['clivi'][:]
    clivi_lat_lon = h5f['clivi_lat_lon'][:]

    clt_lc = h5f['clt_lc'][:]
    clwvi_lc = h5f['clwvi_lc'][:]
    clt_lc_lat_lon = h5f['clt_lc_lat_lon'][:]
    clwvi_lc_lat_lon = h5f['clwvi_lc_lat_lon'][:]

    rsdt_lat_lon = h5f['rsdt_lat_lon'][:]
    rsut_lat_lon = h5f['rsut_lat_lon'][:]
    rsutcs_lat_lon = h5f['rsutcs_lat_lon'][:]
    rtmt_lat_lon = h5f['rtmt_lat_lon'][:]
    albedo_reg = h5f['albedo_reg'][:]
    albedo_so = h5f['albedo_so']

    if model == 'CMIP6-MRI-ESM2':
        rlut_lat_lon = h5f['rlut_lat_lon'][:]
        rlutcs_lat_lon = h5f['rlutcs_lat_lon'][:]

    if model == 'CMIP6-GFDL-AM4':
        mmrdust_lat_lon = h5f['mmrdust_lat_lon'][:]
        mmroa_lat_lon = h5f['mmroa_lat_lon'][:]
        mmrso4_lat_lon = h5f['mmrso4_lat_lon'][:]
        aerosol_norm_lat_lon = h5f['aerosol_norm_lat_lon'][:]


for index,key in enumerate(h5f.keys()):
    print (index, key)
    
#--- sample plots for confirmation ---#

fig, ax = plt.subplots()
ax.plot( constants.lat, clivi )
ax.set_ylabel('Cloud Fraction')
ax.set_xlabel('Latitude')
ax.set_title ('Global Cloud Fraction vs Latitude')
plt.grid(True)
plt.show()


# fig, ax = plt.subplots()
# ax.plot( constants.ta_so, clw_t_so )
# ax.set_ylabel('Cloud Liquid Water Fraction')
# ax.set_xlabel('Temperature')
# ax.set_title ('Cloud Liquid Water Fraction vs Temperature')
# plt.grid(True)
# plt.show()


# fig, ax = plt.subplots()
# ax.plot( cl_g, constants.alt )
# ax.set_ylabel('Altitude (km)')
# ax.set_xlabel('Cloud Liquid Water Fraction')
# ax.set_title ('Southern Ocean Cloud Liquid Water Fraction vs Altitude')
# plt.grid(True)
# plt.show()


# fig, ax = plt.subplots()
# cont = ax.contourf( constants.lat, constants.liq_alt, clw_alt_lat )
# temp = ax.contour( constants.lat, constants.liq_alt, (ta_alt_lat - 273.15), colors='white')
# temp.collections[5].set_linewidth(3)
# temp.collections[5].set_color('white')
# ax.clabel(temp, inline=1, fontsize=10)
# ax.set_xlabel('Latitude')
# ax.set_ylabel('Altitude (km)')
# cbar = fig.colorbar(cont, orientation='horizontal')
# cbar.set_label('Cloud liquid Water Fraction')
# plt.show()

if model == 'CERES' or model == 'CALIPSO':
    lon = constants.lon - 180
else:
    lon = constants.lon
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
ax.coastlines()
p = ax.contourf(lon, constants.lat, np.transpose(albedo_reg), transform=ccrs.PlateCarree(), cmap='coolwarm')
cbar = plt.colorbar(p, orientation='horizontal')
cbar.set_label('Cloud Fraction')
ax.set_title('Total Cloud Fraction - CMIP6-GFDL-AM4')

plt.show()




