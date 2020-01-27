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
# import cartopy.crs as ccrs

#--- Set Location, date period and model ---#

# specify model from the list above
model = 'CALIPSO' 

# specify location: home, uni, hdd, laptop
location = constants.home + 'climate-analysis/reduced_data'

os.chdir( location )

if model == 'CCCM':
    h5f = h5py.File( constants.date_cccm + '_' + model + '.h5', 'r')
    cl_alt_lat = h5f['cl_alt_lat'][:]
    cl_g = h5f['cl_g'][:]
    cl_so = h5f['cl_so'][:]
    clic_alt_lat = h5f['cli_alt_lat'][:]
    cli_g = h5f['cli_g'][:]
    cli_so = h5f['cli_so'][:]
    clwc_alt_lat = h5f['clw_alt_lat'][:]
    clw_frac_alt_lat = h5f['clw_frac_alt_lat'][:]
    clw_g = h5f['clw_g'][:]
    clw_so = h5f['clw_so'][:]
    ta_alt_lat = h5f['ta_alt_lat'][:]
    clwc_t_g = h5f['clw_t_g'][:]
    clwc_t_so = h5f['clw_t_so'][:]
    clw_frac_t_g = h5f['clw_frac_t_g'][:]
    clw_frac_t_so = h5f['clw_frac_t_so'][:]
    cl_t_g = h5f['cl_t_g'][:]
    cl_t_so = h5f['cl_t_so'][:]
    full_clw_alt_lat = h5f['full_clw_alt_lat'][:]
    full_ta_alt_lat = h5f['full_ta_alt_lat'][:]
    # clt = h5f['clt'][:]
    # clt_lat_lon = h5f['clt_lat_lon'][:]
    # clwvi = h5f['clwvi'][:]
    # clwvi_lat_lon = h5f['clwvi_lat_lon'][:]
    # clivi = h5f['clivi'][:]
    # clivi_lat_lon = h5f['clivi_lat_lon'][:]
   
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
 
    w_alt_lat = h5f['w_alt_lat'][:]

if model == 'CALIPSO':    
    h5f = h5py.File( constants.date_cmip6 + '_' + 'CALIPSO.h5', 'r')
    cl_alt_lat = h5f['cl_alt_lat'][:]
    cl_g = h5f['cl_g'][:]
    cl_so = h5f['cl_so'][:]
    cli_frac_alt_lat = h5f['cli_frac_alt_lat'][:]
    cli_frac_g = h5f['cli_frac_g'][:]
    cli_frac_so = h5f['cli_frac_so'][:]
    clw_frac_g = h5f['clw_frac_g'][:]
    clw_frac_so = h5f['clw_frac_so'][:]
    clw_frac_t_g = h5f['clw_frac_t_g'][:]
    clw_frac_t_so = h5f['clw_frac_t_so'][:]
    clw_frac_alt_lat = h5f['clw_frac_alt_lat'][:]
    cl_t_g = h5f['cl_t_g'][:]
    cl_t_so = h5f['cl_t_so'][:]
    clt = h5f['clt'][:]
    clt_lat_lon = h5f['clt_lat_lon'][:]
    clw_frac = h5f['clw_frac'][:]
    clw_frac_lat_lon = h5f['clw_frac_lat_lon'][:]
    cli_frac = h5f['cli_frac'][:]
    cli_frac_lat_lon = h5f['cli_frac_lat_lon'][:]
    full_clw_frac_alt_lat = h5f['full_clw_frac_alt_lat'][:]
    full_ta_alt_lat = h5f['full_ta_alt_lat'][:]
    ta_alt_lat = h5f['ta_alt_lat'][:]
    clt_l = h5f['clt_l'][:]
    clw_frac_l = h5f['clw_frac_l'][:]
    clt_l_lat_lon = h5f['clt_l_lat_lon'][:]
    clw_frac_l_lat_lon = h5f['clw_frac_l_lat_lon'][:]

    
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
    clwc_alt_lat = h5f['clwc_alt_lat'][:]
    clw_g = h5f['clw_g'][:]
    clw_so = h5f['clw_so'][:]
    clwc_g = h5f['clwc_g'][:]
    clwc_so = h5f['clwc_so'][:]
    ta_alt_lat = h5f['ta_alt_lat'][:]
    clw_t_g = h5f['clw_t_g'][:]
    clw_t_so = h5f['clw_t_so'][:]
    clwc_t_g = h5f['clwc_t_g'][:]
    clwc_t_so = h5f['clwc_t_so'][:]
    full_clw_alt_lat = h5f['full_clw_alt_lat'][:]
    full_ta_alt_lat = h5f['full_ta_alt_lat'][:]
    clt = h5f['clt'][:]
    clt_lat_lon = h5f['clt_lat_lon'][:]
    clwvi = h5f['clwvi'][:]
    clwvi_lat_lon = h5f['clwvi_lat_lon'][:]
    clivi = h5f['clivi'][:]
    clivi_lat_lon = h5f['clivi_lat_lon'][:]
    full_air_density_alt_lat = h5f['full_air_density_alt_lat'][:]
    liq_air_density_alt_lat = h5f['liq_air_density_alt_lat'][:]

    clt_lc = h5f['clt_lc'][:]
    clt_lc_lat_lon = h5f['clt_lc_lat_lon'][:]
    clwvi_lc_lat_lon = h5f['clwvi_lc_lat_lon'][:]

    rtmt_lat_lon = h5f['rtmt_lat_lon'][:]
    albedo_reg = h5f['albedo_reg'][:]

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
    clw_frac_alt_lat = h5f['clw_frac_alt_lat'][:]
    clwc_alt_lat = h5f['clwc_alt_lat'][:]
    clw_g = h5f['clw_g'][:]
    clw_so = h5f['clw_so'][:]
    clwc_g = h5f['clwc_g'][:]
    clwc_so = h5f['clwc_so'][:]
    ta_alt_lat = h5f['ta_alt_lat'][:]
    clw_t_g = h5f['clw_t_g'][:]
    clw_t_so = h5f['clw_t_so'][:]
    clw_frac_t_g = h5f['clw_frac_t_g'][:]
    clw_frac_t_so = h5f['clw_frac_t_so'][:]
    clwc_t_g = h5f['clwc_t_g'][:]
    clwc_t_so = h5f['clwc_t_so'][:]
    full_clw_alt_lat = h5f['full_clw_alt_lat'][:]
    full_ta_alt_lat = h5f['full_ta_alt_lat'][:]
    clt = h5f['clt'][:]
    clt_lat_lon = h5f['clt_lat_lon'][:]
    clwvi = h5f['clwvi'][:]
    clwvi_lat_lon = h5f['clwvi_lat_lon'][:]
    clivi = h5f['clivi'][:]
    clivi_lat_lon = h5f['clivi_lat_lon'][:]
    full_air_density_alt_lat = h5f['full_air_density_alt_lat'][:]
    liq_air_density_alt_lat = h5f['liq_air_density_alt_lat'][:]

    clt_lc = h5f['clt_lc'][:]
    clt_lc_lat_lon = h5f['clt_lc_lat_lon'][:]
    clwvi_lc_lat_lon = h5f['clwvi_lc_lat_lon'][:]

    rtmt_lat_lon = h5f['rtmt_lat_lon'][:]
    albedo_reg = h5f['albedo_reg'][:]

    if model == 'CMIP6-GFDL-AM4':
        mmrdust_lat_lon = h5f['mmrdust_lat_lon'][:]
        mmroa_lat_lon = h5f['mmroa_lat_lon'][:]
        mmrso4_lat_lon = h5f['mmrso4_lat_lon'][:]
        aerosol_norm_lat_lon = h5f['aerosol_norm_lat_lon'][:]


for index,key in enumerate(h5f.keys()):
    print (index, key)

# cl_alt_lat = np.transpose(cl_alt_lat)
# clwc_alt_lat = np.transpose(clwc_alt_lat)
# lat = constants.lat
# print(constants.globalalt_latMeanVal(cl_alt_lat, constants.lat))
# print(constants.globalalt_latMeanVal(cl_alt_lat[constants.so_idx_1:constants.so_idx_2], constants.lat[constants.so_idx_1:constants.so_idx_2]))
# print(constants.globalalt_latMeanVal(clwc_alt_lat, constants.lat)

# print(constants.globalalt_latMeanVal(clwc_alt_lat[constants.so_idx_1:constants.so_idx_2], lat[constants.so_idx_1:constants.so_idx_2]))


#--- sample plots for confirmation ---#

fig, ax = plt.subplots()
ax.plot( constants.lat[constants.lat_confine_1:constants.lat_confine_2], clt[constants.lat_confine_1:constants.lat_confine_2] )
ax.set_ylabel('Cloud Fraction')
ax.set_xlabel('Latitude')
ax.set_title ('Global Cloud Fraction vs Latitude')
plt.grid(True)
plt.show()

# fig, ax = plt.subplots()
# ax.plot( constants.lat[constants.lat_confine_1:constants.lat_confine_2], clwvi[constants.lat_confine_1:constants.lat_confine_2] )
# ax.plot( constants.lat[constants.lat_confine_1:constants.lat_confine_2], clivi[constants.lat_confine_1:constants.lat_confine_2] )
# ax.set_ylabel('Cloud Fraction')
# ax.set_xlabel('Latitude')
# ax.set_title ('Global Cloud Fraction vs Latitude')
# plt.grid(True)
# plt.show()

# fig, ax = plt.subplots()
# ax.plot( constants.ta, cl_t_g )
# ax.plot( constants.ta, cl_t_so )
# ax.set_ylabel('Mean Cloud Liquid Water Mass Fraction in Air (kg/kg)')
# ax.set_xlabel('Temperature (K)')
# ax.set_title ('Mean Cloud Liquid Water Mass Fraction vs Temperature')
# ax.axvline(x=273, label = '273K', color = 'black', linestyle='--')

# plt.grid(True)
# plt.show()


# fig, ax = plt.subplots()
# ax.plot( cl_so, constants.alt )
# ax.plot( cl_g, constants.alt )
# ax.set_ylabel('Altitude (km)')
# ax.set_xlabel('Mean Cloud  Fraction ')
# ax.set_title ('Cloud Fraction vs Altitude')
# plt.grid(True)
# plt.show()

# fig, ax = plt.subplots()
# ax.plot( clwc_so, constants.liq_alt )
# ax.plot( clwc_g, constants.liq_alt )
# ax.plot( cli_so, constants.alt )
# ax.plot( cli_g, constants.alt )
# ax.set_ylabel('Altitude (km)')
# ax.set_xlabel('Mean Cloud Liquid Water Mass Fraction in Air (kg/kg)')
# ax.set_title ('Southern Ocean Cloud Liquid Water Fraction vs Altitude')
# plt.grid(True)
# plt.show()


# fig, ax = plt.subplots()
# cont = ax.contourf( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.alt, full_clw_frac_alt_lat[:,constants.lat_confine_1:constants.lat_confine_2] )
# temp = ax.contour( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.alt, (full_ta_alt_lat[:,constants.lat_confine_1:constants.lat_confine_2] - 273.15), colors='white')
# temp.collections[5].set_linewidth(3)
# temp.collections[5].set_color('white')
# ax.clabel(temp, inline=1, fontsize=10)
# ax.set_xlabel('Latitude')
# ax.set_ylabel('Altitude (km)')
# cbar = fig.colorbar(cont, orientation='horizontal')
# cbar.set_label('Mean Cloud Liquid Water Mass Fraction in Air (kg/kg)')
# plt.show()


# ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
# ax.coastlines()
# p = ax.contourf(constants.lon, constants.lat, clt_lat_lon, transform=ccrs.PlateCarree(), cmap='coolwarm')
# cbar = plt.colorbar(p, orientation='horizontal')
# cbar.set_label('Cloud Fraction')
# ax.set_title('Total Cloud Fraction - ' + model)
# plt.show()