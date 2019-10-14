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
    MISR
"""

#read hdf5
import h5py
import os
import matplotlib.pyplot as plt
import constants
import numpy as np
import cartopy.crs as ccrs

#--- Set Location, date period and model ---#

model = 'CMIP5-GFDL-HIRAM-C360'
location = constants.home + 'climate-analysis/reduced_data'
date_cmip5 = 'Jan_2001_Dec_2005'
date_cmip6 = 'Jan_2006_Dec_2010'
date_ceres = 'Jun_2006_Jun_2011'
os.chdir( location )

if model == 'CCCM':
    h5f = h5py.File( 'Jun_2006_Apr_2011_CCCM.h5', 'r')
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
   
elif model == 'CERES':
    h5f = h5py.File( date_ceres + '_' + model + '.h5', 'r')
    clt = h5f['clt'][:]
    clt_lat_lon = h5f['clt_lat_lon'][:]
    clwvi = h5f['clwvi'][:]
    clwvi_lat_lon = h5f['clwvi_lat_lon'][:]
    clivi = h5f['clivi'][:]
    clivi_lat_lon = h5f['clivi_lat_lon'][:]
#    clr_toa_sw_reg = h5f['clr_toa_sw_reg'][:]
#    clr_toa_sw_reg = h5f['clr_toa_sw_reg'][:]
#    clr_toa_lw_reg = h5f['clr_toa_lw_reg'][:]
#    clr_toa_lw_reg = h5f['clr_toa_lw_reg'][:]
#    clr_toa_net_reg = h5f['clr_toa_net_reg'][:]
#    clr_toa_net_reg = h5f['clr_toa_net_reg'][:]
#    all_toa_sw_reg = h5f['all_toa_sw_reg'][:]
#    all_toa_sw_reg = h5f['all_toa_sw_reg'][:]
#    all_toa_lw_reg = h5f['all_toa_lw_reg'][:]
#    all_toa_lw_reg_so = h5f['all_toa_lw_reg_so'][:]
#    all_toa_net_reg = h5f['all_toa_net_reg'][:]
#    all_toa_net_reg_so = h5f['all_toa_net_reg_so'][:]
    

elif model == 'ECMWF':    
    h5f = h5py.File( date_cmip6 + '_' + model + '.h5', 'r')
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

elif model == 'CALIPSO':    
    h5f = h5py.File( 'Jun_2006_Jun_2011_CALIPSO.h5', 'r')
    cal_cl_alt_lat = h5f['cl_alt_lat'][:]
    cal_cl_g = h5f['cl_g'][:]
    cal_cl_so = h5f['cl_so'][:]
    cal_cli_alt_lat = h5f['cli_alt_lat'][:]
    cal_cli_g = h5f['cli_g'][:]
    cal_cli_so = h5f['cli_so'][:]
    cal_clw_g = h5f['clw_g'][:]
    cal_clw_so = h5f['clw_so'][:]
    cal_clw_t_g = h5f['clw_t_g'][:]
    cal_clw_t_so = h5f['clw_t_so'][:]
    cal_clw_alt_lat = h5f['clw_alt_lat'][:]
    cal_clt = h5f['clt'][:]
    cal_clt_lat_lon = h5f['clt_lat_lon'][:]
    cal_clwvi = h5f['clwvi'][:]
    cal_clwvi_lat_lon = h5f['clwvi_lat_lon'][:]
    cal_clivi = h5f['clivi'][:]
    cal_clivi_lat_lon = h5f['clivi_lat_lon'][:]
    cal_full_clw_alt_lat = h5f['full_clw_alt_lat'][:]
#    cal_full_ta_alt_lat = h5f['full_ta_alt_lat'][:]
#    cal_ta_alt_lat = h5f['ta_alt_lat'][:]

elif model == 'MISR':    
    h5f = h5py.File( 'Jan_2006_Dec_2010_MISR.h5', 'r')
    cl_g = h5f['cl_g'][:]
    cl_so = h5f['cl_so'][:]
    
elif model == 'CMIP5-CESM1-CAM5' or model == 'CMIP5-GFDL-HIRAM-C360' or model == 'CMIP5-GISS-E2R' or model == 'CMIP5-IPSL-CM5A-LR' or model == 'CMIP5-MIROC5' or model == 'CMIP5-MRI-CGCM3': 
    h5f = h5py.File( date_cmip5 + '_' + model + '.h5', 'r')
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


else:
    h5f = h5py.File( date_cmip6 + '_' + model + '.h5', 'r')
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


fig, ax = plt.subplots()
ax.plot( constants.ta_so, clw_t_so )
ax.set_ylabel('Cloud Liquid Water Fraction')
ax.set_xlabel('Temperature')
ax.set_title ('Cloud Liquid Water Fraction vs Temperature')
plt.grid(True)
plt.show()


fig, ax = plt.subplots()
ax.plot( clw_so, constants.liq_alt )
ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Cloud Liquid Water Fraction')
ax.set_title ('Southern Ocean Cloud Liquid Water Fraction vs Altitude')
plt.grid(True)
plt.show()


fig, ax = plt.subplots()
cont = ax.contourf( constants.lat, constants.liq_alt, clw_alt_lat )
temp = ax.contour( constants.lat, constants.liq_alt, (ta_alt_lat - 273.15), colors='white')
temp.collections[5].set_linewidth(3)
temp.collections[5].set_color('white')
ax.clabel(temp, inline=1, fontsize=10)
ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cont, orientation='horizontal')
cbar.set_label('Cloud liquid Water Fraction')
plt.show()


ax = plt.axes(projection=ccrs.Mollweide(central_longitude=180))
ax.coastlines()
p = ax.contourf(constants.lon, constants.lat, clt_lat_lon, transform=ccrs.PlateCarree())
cbar = plt.colorbar(p, orientation='horizontal')
cbar.set_label('Cloud Fraction')
plt.show()




