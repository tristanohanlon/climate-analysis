# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 18:40:37 2019

@author: Tristan O'Hanlon - University of Auckland & Jonathan Rogers

These functions will create a standard model dataset when parsed:
    
reduce_dataset( directory, filename, save_location, start, end, cl_is_fractional=False )

example:
reduce_dataset( 'gfdl_am4', '_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'E:/University/University/MSc/Models/climate-analysis/GFDL-AM4', datetime.datetime( 2000, 1, 1 ), datetime.datetime( 2006, 1, 1 ) )

use ECMWF ta_alt_lat over contour data
Specify time range at the bottom

Change lon in constants from 0-360 to -180 to 180
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
from scipy import ndimage as nd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


location = constants.home # home, uni, hdd or laptop
os.chdir( location + 'Data/CALIPSO-GOCCP' )

f = Dataset('MapLowMidHigh330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')
raw_lat = f.variables['latitude'][4:85]
start_idx = np.abs(raw_lat - (-70)).argmin()
end_idx = np.abs(raw_lat - (-50)).argmin()

##### IMPORTANT - change lon in constants to -180 to 180
raw_lon = f.variables['longitude'][:] - 180

#opaque (thick) clouds
cc = np.array(f.variables['cltcalipso'][:48]) # 7.2006 to 12.2020 - 54 months
cc[cc < 0] = None #set fill values to nan
clt_lat_lon  = np.roll(np.nanmean(cc, axis = 0)[4:85], 179) #average over time


print(constants.global2DMean(clt_lat_lon[start_idx:end_idx], raw_lat[start_idx:end_idx]))




# clt = np.nanmean(clt_lat_lon, axis = -1) #average over longitude (lat)


# ccl = np.array(f.variables['cllcalipso'][5:55]) # 7.2006 to 12.2020 - 54 months
# ccl[ccl < 0] = None #set fill values to nan
# clt_lc_lat_lon  = np.roll(np.nanmean(ccl, axis = 0)[4:85], 179) #average over time
# clt_lc_lat_lon = constants.fill(clt_lc_lat_lon)
# clt_lc = np.nanmean(clt_lc_lat_lon, axis = -1) #average over longitude (lat)



# ###########################---get lat-lon - phase fraction---###########################

# f = Dataset('MapLowMidHigh_Phase330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')
# clwvi = np.array(f.variables['cltcalipso_liq'][5:55]) # 7.2006 to 12.2020 - 54 months
# clwvi[clwvi < 0] = None #set fill values to nan
# clwvi_lat_lon = np.roll(np.nanmean(clwvi, axis = 0)[4:85], 179) #average over time
# clwvi = np.nanmean(clwvi_lat_lon, axis = -1) #average over longitude

# clivi = np.array(f.variables['cltcalipso_ice'][5:55]) # 7.2006 to 12.2020 - 54 months
# clivi[clivi < 0] = None #set fill values to nan
# clivi_lat_lon = np.roll(np.nanmean(clivi, axis = 0)[4:85], 179) #average over time
# clivi = np.nanmean(clivi_lat_lon, axis = -1) #average over longitude

# clwvi_lc = np.array(f.variables['cllcalipso_liq'][5:55]) # 7.2006 to 12.2020 - 54 months
# clwvi_lc[clwvi_lc < 0] = None #set fill values to nan
# clwvi_lc_lat_lon = np.roll(np.nanmean(clwvi_lc, axis = 0)[4:85], 179) #average over time
# clwvi_lc = np.nanmean(clwvi_lc_lat_lon, axis = -1) #average over longitude

# clivi_lc = np.array(f.variables['cllcalipso_ice'][5:55]) # 7.2006 to 12.2020 - 54 months
# clivi_lc[clivi_lc < 0] = None #set fill values to nan
# clivi_lc_lat_lon = np.roll(np.nanmean(clivi_lc, axis = 0)[4:85], 179) #average over time
# clivi_lc = np.nanmean(clivi_lc_lat_lon, axis = -1) #average over longitude


# ###########################---get alt - cloud fraction---###########################

# f = Dataset('3D_CloudFraction330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')
# raw_alt = np.array(f.variables['alt_mid'][:]) # km

# cl = np.array(f.variables['clcalipso'][5:55]) # 7.2006 to 12.2020 - 54 months
# cl[cl < 0] = None #set fill values to nan
# cl = constants.fill(np.nanmean(cl, axis = 0)[:,4:85]) #average over time
# cl_alt_lat = np.nanmean(cl, axis = -1) #average over longitude

# cl_g = constants.global3DMean(cl, raw_lat)
# cl_so = constants.global3DMean(cl[:,start_idx:end_idx], raw_lat[start_idx:end_idx])

# #ax = plt.axes(projection=ccrs.Mollweide(central_longitude=180))
# #ax.coastlines()
# #p = ax.contourf(raw_lon, raw_lat, clw_lc_lat_lon, transform=ccrs.PlateCarree(), cmap='coolwarm')
# #cbar = plt.colorbar(p, orientation='horizontal')
# #cbar.set_label('Cloud Fraction')
# #ax.set_title('Total Cloud Fraction - CMIP6-GFDL-AM4')
# #
# #plt.show()

# ###########################---get alt - phase fractions---###########################

# f = Dataset('3D_CloudFraction_Phase330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')

# #liquid water fraction
# # clw = np.array(f.variables['clcalipso_liq'][5:55]) # 5 years
# # clw[clw < 0] = None #set fill values to nan
# # clw = constants.fill(np.nanmean(clw, axis = 0)[:,4:85]) #average over time
# # clw_alt_lat = np.nanmean(clw, axis = -1) #average over longitude

# #liquid water fraction
# cli = np.array(f.variables['clcalipso_RPIC'][5:55]) # 5 years
# cli[cli < 0] = None #set fill values to nan
# cli = constants.fill(np.nanmean(cli, axis = 0)[:,4:85]) #average over time
# cli_alt_lat = np.nanmean(cli, axis = -1) #average over longitude

# clw = 1 - cli
# clw_alt_lat = (1 - cli_alt_lat) * cl_alt_lat
# cli_alt_lat = cli_alt_lat * cl_alt_lat 

# clw_g = constants.global3DMean(clw, raw_lat)
# clw_so = constants.global3DMean(clw[:,start_idx:end_idx], raw_lat[start_idx:end_idx])
# cli_g = constants.global3DMean(cli, raw_lat)
# cli_so = constants.global3DMean(cli[:,start_idx:end_idx], raw_lat[start_idx:end_idx])


# ###########################---get alt - temp - phase fractions---###########################

# f = Dataset('3D_CloudFraction_Temp330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')
# raw_ta = np.array(f.variables['temp_mid'][:]) + 273 # km (38)

# cl_t = np.array(f.variables['cltemp'][5:55]) # 7.2006 to 12.2020 - 54 months
# cl_t [cl_t  < 0] = None #set fill values to nan
# cl_t  = constants.fill(np.nanmean(cl_t, axis = 0)[:,4:85]) #average over time
# cl_t_lat =  np.nanmean(cl_t, axis = -1) #average over longitude

# clw_t = np.array(f.variables['cltemp_liq'][5:55]) # 7.2006 to 12.2020 - 54 months
# clw_t[clw_t < 0] = None #set fill values to nan
# clw_t = constants.fill(np.nanmean(clw_t, axis = 0)[:,4:85]) #average over time
# clw_t_lat =  np.nanmean(clw_t, axis = -1) #average over longitude

# cl_t_g = constants.global3DMean(cl_t, raw_lat)
# cl_t_so = constants.global3DMean(cl_t[:,start_idx:end_idx], raw_lat[start_idx:end_idx])
# clw_t_g = constants.global3DMean(clw_t, raw_lat)
# clw_t_so = constants.global3DMean(clw_t[:,start_idx:end_idx], raw_lat[start_idx:end_idx])


# ######################################

# interpolated = interpolate.interp1d(raw_lat, clt, kind = 'cubic', fill_value = 'extrapolate')
# clt = interpolated(constants.lat)

# interpolated = interpolate.interp1d(raw_lat, clt_lc, kind = 'cubic', fill_value = 'extrapolate')
# clt_lc = interpolated(constants.lat)

# interpolated = interpolate.interp1d(raw_lat, clwvi, kind = 'cubic', fill_value = 'extrapolate')
# clwvi = interpolated(constants.lat)

# interpolated = interpolate.interp1d(raw_lat, clwvi_lc, kind = 'cubic', fill_value = 'extrapolate')
# clwvi_lc = interpolated(constants.lat)

# interpolated = interpolate.interp1d(raw_lat, clivi, kind = 'cubic', fill_value = 'extrapolate')
# clivi = interpolated(constants.lat)

# interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clt_lat_lon ), fill_value = np.nan)
# clt_lat_lon = interpolated(constants.lat, constants.lon)
# clt_lat_lon = np.transpose(clt_lat_lon)

# interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clwvi_lat_lon ), kind = 'cubic', fill_value = np.nan)
# clwvi_lat_lon = interpolated(constants.lat, constants.lon)
# clwvi_lat_lon = np.transpose(clwvi_lat_lon)

# interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clt_lc_lat_lon ), fill_value = np.nan)
# clt_lc_lat_lon = interpolated(constants.lat, constants.lon)
# clt_lc_lat_lon = np.transpose(clt_lc_lat_lon)

# interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clwvi_lc_lat_lon ), kind = 'cubic', fill_value = np.nan)
# clwvi_lc_lat_lon = interpolated(constants.lat, constants.lon)
# clwvi_lc_lat_lon = np.transpose(clwvi_lc_lat_lon)

# interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clivi_lat_lon ), kind = 'cubic', fill_value = np.nan)
# clivi_lat_lon = interpolated(constants.lat, constants.lon)
# clivi_lat_lon = np.transpose(clivi_lat_lon)

# interpolated = interpolate.interp1d(raw_alt, cl_g, kind = 'cubic', fill_value = 'extrapolate')
# cl_g = interpolated(constants.alt)

# interpolated = interpolate.interp1d(raw_alt, cl_so, kind = 'cubic', fill_value = 'extrapolate')
# cl_so = interpolated(constants.alt)

# interpolated = interpolate.interp1d(raw_alt, clw_g, kind = 'cubic', fill_value = 'extrapolate')
# clw_g = interpolated(constants.liq_alt)

# interpolated = interpolate.interp1d(raw_alt, cli_g, kind = 'cubic', fill_value = 'extrapolate')
# cli_g = interpolated(constants.alt)

# interpolated = interpolate.interp1d(raw_alt, clw_so, kind = 'cubic', fill_value = 'extrapolate')
# clw_so = interpolated(constants.liq_alt)

# interpolated = interpolate.interp1d(raw_alt, cli_so, kind = 'cubic', fill_value = 'extrapolate')
# cli_so = interpolated(constants.alt)

# interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( cl_alt_lat ), kind = 'cubic', fill_value = np.nan)
# cl_alt_lat = interpolated(constants.alt, constants.lat)

# interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( clw_alt_lat ), kind = 'cubic', fill_value = np.nan)
# full_clw_alt_lat = interpolated(constants.alt, constants.lat)

# interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( clw_alt_lat ), kind = 'cubic', fill_value = np.nan)
# clw_alt_lat = interpolated(constants.liq_alt, constants.lat)

# interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( cli_alt_lat ), kind = 'cubic', fill_value = np.nan)
# cli_alt_lat = interpolated(constants.alt, constants.lat)
# cli_alt_lat [cli_alt_lat  < 0] = None #set fill values to nan

# interpolated = interpolate.interp1d(raw_ta, cl_t_g, kind = 'cubic', fill_value = np.nan)
# cl_t_g = interpolated(constants.ta)

# interpolated = interpolate.interp1d(raw_ta, cl_t_so, kind = 'cubic', fill_value = np.nan)
# cl_t_so = interpolated(constants.ta)

# interpolated = interpolate.interp1d(raw_ta, clw_t_g, kind = 'cubic', fill_value = np.nan)
# clw_t_g = interpolated(constants.ta)

# interpolated = interpolate.interp1d(raw_ta, clw_t_so, kind = 'cubic', fill_value = np.nan)
# clw_t_so = interpolated(constants.ta)


# fig, ax = plt.subplots()
# ax.plot( constants.ta, clw_t_g )
# ax.plot( constants.ta, clw_t_so )
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
# ax.set_xlabel('Mean Cloud Fraction ')
# ax.set_title ('Cloud Fraction vs Altitude')
# plt.grid(True)
# plt.show()

# fig, ax = plt.subplots()
# ax.plot( clw_so, constants.liq_alt )
# ax.plot( clw_g, constants.liq_alt )
# ax.plot( cli_so, constants.alt )
# ax.plot( cli_g, constants.alt )
# ax.set_ylabel('Altitude (km)')
# ax.set_xlabel('Mean Cloud Liquid Water Mass Fraction in Air (kg/kg)')
# ax.set_title ('Southern Ocean Cloud Liquid Water Fraction vs Altitude')
# plt.grid(True)
# plt.show()


# fig, ax = plt.subplots()
# cont = ax.contourf( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.liq_alt, np.transpose(clw_alt_lat[constants.lat_confine_1:constants.lat_confine_2]) )
# ax.set_xlabel('Latitude')
# ax.set_ylabel('Altitude (km)')
# cbar = fig.colorbar(cont, orientation='horizontal')
# cbar.set_label('Mean Cloud Liquid Water Mass Fraction in Air (kg/kg)')
# plt.show()



#     #----------------------------#

# os.chdir(location + 'climate-analysis/reduced_data')

# with h5py.File( 'Jan_2007_Dec_2010_CALIPSO.h5', 'w' ) as p:
    
#     p.create_dataset('clt', data=clt) # total cloud fraction corresponding to lat
#     p.create_dataset('clt_lat_lon', data=clt_lat_lon ) # total cloud fraction corresponding to lat, lon
  
#     p.create_dataset('clt_lc', data=clt_lc) # total cloud fraction corresponding to lat
#     p.create_dataset('clt_lc_lat_lon', data=clt_lc_lat_lon ) # total cloud fraction corresponding to lat, lon    
    
#     p.create_dataset('clwvi', data=clwvi) # total cloud liquid water fraction corresponding to lat
#     p.create_dataset('clwvi_lat_lon', data=clwvi_lat_lon ) # total cloud fraction corresponding to lat, lon

#     p.create_dataset('clwvi_lc', data=clwvi_lc) # total cloud fraction corresponding to lat
#     p.create_dataset('clwvi_lc_lat_lon', data=clwvi_lc_lat_lon ) # total cloud fraction corresponding to lat, lon    

#     p.create_dataset('clivi', data=clivi) # total cloud ice fraction corresponding to lat
#     p.create_dataset('clivi_lat_lon', data=clivi_lat_lon ) # total cloud fraction corresponding to lat, lon
    
#     p.create_dataset('cl_g', data=cl_g) # global layer total cloud fraction corresponding to alt
#     p.create_dataset('clw_frac_g', data=clw_g) # global layer cloud liquid water fraction corresponding to liq_alt
#     p.create_dataset('cli_frac_g', data=cli_g) # global layer cloud ice water fraction corresponding to alt
    
#     p.create_dataset('cl_so', data=cl_so) # southern ocean layer total cloud fraction corresponding to alt
#     p.create_dataset('clw_frac_so', data=clw_so) # southern ocean layer cloud liquid water fraction corresponding to liq_alt
#     p.create_dataset('cli_frac_so', data=cli_so) # southern ocean layer cloud ice water fraction corresponding to alt

#     p.create_dataset('cl_t_g', data=cl_t_g) # global layer cloud liquid water fraction corresponding to ta_g
#     p.create_dataset('cl_t_so', data=cl_t_so) # global layer cloud liquid water fraction corresponding to ta_so

#     p.create_dataset('clw_frac_t_g', data=clw_t_g) # global layer cloud liquid water fraction corresponding to ta_g
#     p.create_dataset('clw_frac_t_so', data=clw_t_so) # global layer cloud liquid water fraction corresponding to ta_so
   
#     p.create_dataset('cl_alt_lat', data=np.transpose(cl_alt_lat) ) # total cloud fraction corresponding to alt and lat
#     p.create_dataset('cli_frac_alt_lat', data=np.transpose(cli_alt_lat) ) # cloud ice water fraction corresponding to alt and lat
#     p.create_dataset('clw_frac_alt_lat', data=np.transpose(clw_alt_lat) )  # cloud liquid water fraction corresponding to liq_alt and lat
#     p.create_dataset('full_clw_frac_alt_lat', data=np.transpose(full_clw_alt_lat) )  # cloud liquid water fraction corresponding to liq_alt and lat


#     p.close()


#     take from ECMWF - import first before saving
    
"""    
import numpy as np
import os
import h5py
import constants

location = constants.home + 'climate-analysis/reduced_data'
os.chdir( location )
h5f = h5py.File( 'Jan_2006_Dec_2010_ECMWF.h5', 'r')
ta_alt_lat = h5f['ta_alt_lat'][:]
full_ta_alt_lat = h5f['full_ta_alt_lat'][:]

f = h5py.File('Jan_2007_Dec_2010_CALIPSO.h5','a')
f.create_dataset('full_ta_alt_lat', data=full_ta_alt_lat)
f.create_dataset('ta_alt_lat', data=ta_alt_lat)
f.close()
"""
