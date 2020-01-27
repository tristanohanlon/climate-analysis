# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon - University of Auckland & Jonathan Rogers


"""
import numpy as np
import os
from netCDF4 import Dataset
from netCDF4 import date2index
import h5py
import math
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline
from scipy import stats
import datetime
from pprint import pprint
from sklearn.impute import SimpleImputer
import constants
import time
from scipy import ndimage as nd
import matplotlib.pyplot as plt

###############################################################################
start = time.time()

#set location
location = constants.home
os.chdir( location + 'Data/ECMWF/' )

dataset = Dataset('200601-201012_ECMWF_amon_clt_cltlc_clwvi_clivi.nc', 'r')

#lon = dataset.variables['longitude'][:] #Extract longitude data
raw_lat = np.flip(np.array(dataset.variables['latitude'][:]), axis = 0) #Extract latitude data
raw_lon = np.array(dataset.variables['longitude'][:]) #Extract latitude data

clt = np.array(dataset.variables['tcc'][:]) #Extract total cloud cover, keyed to time, lon and lat
clt_lat_lon = np.flip(np.mean(clt, axis = 0), axis = 0)

clt_lc = np.array(dataset.variables['lcc'][:]) #Extract low cloud cover, keyed to time, lon and lat
clt_lc_lat_lon = np.flip(np.mean(clt_lc, axis = 0), axis = 0)

clivi = np.array(dataset.variables['tciw'][:]) #Extract Ice water path (kg/m^2), keyed to time, lon and lat
clivi_lat_lon = np.flip(np.mean(clivi, axis = 0), axis = 0)

clwvi = np.array(dataset.variables['tclw'][:]) #Extract liquid water path (kg/m^2), keyed to time, lon and lat
clwvi_lat_lon = np.flip(np.mean(clwvi, axis = 0), axis = 0)

clwvi = np.nanmean(clwvi_lat_lon , axis = -1)
clivi = np.nanmean(clivi_lat_lon , axis = -1)
clt = np.nanmean(clt_lat_lon , axis = -1)
clt_lc = np.nanmean(clt_lc_lat_lon , axis = -1)

interpolated = interpolate.interp1d(raw_lat, clt, kind = 'cubic')
clt = interpolated(constants.lat)

interpolated = interpolate.interp1d(raw_lat, clt_lc, kind = 'cubic')
clt_lc = interpolated(constants.lat)

interpolated = interpolate.interp1d(raw_lat, clwvi, kind = 'cubic')
clwvi = interpolated(constants.lat)

interpolated = interpolate.interp1d(raw_lat, clivi, kind = 'cubic')
clivi = interpolated(constants.lat)

interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clt_lat_lon ), kind = 'cubic')
clt_lat_lon = interpolated(constants.lat, constants.lon)
clt_lat_lon = np.transpose(clt_lat_lon)

interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clt_lc_lat_lon ), kind = 'cubic')
clt_lc_lat_lon = interpolated(constants.lat, constants.lon)
clt_lc_lat_lon = np.transpose(clt_lc_lat_lon)

interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clwvi_lat_lon ), kind = 'cubic')
clwvi_lat_lon = interpolated(constants.lat, constants.lon)
clwvi_lat_lon = np.transpose(clwvi_lat_lon)

interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clivi_lat_lon ), kind = 'cubic')
clivi_lat_lon = interpolated(constants.lat, constants.lon)
clivi_lat_lon = np.transpose(clivi_lat_lon)


#------------- Profile variables -------------#


with Dataset( '200601-201012_ECMWF_plevel_T_cc_clw_cli_w.nc', 'r') as f:
    p = np.flip( constants.extract_data( f, 'level' ), axis = 0 )        
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
    raw_lat = np.flip( constants.extract_data( f, 'latitude' ), axis = 0 )
    start_idx = np.abs(raw_lat - (-70)).argmin()
    end_idx = np.abs(raw_lat - (-50)).argmin()



    cl = np.flip( np.flip( np.mean( constants.extract_data( f, 'cc' ), axis = 0 ), axis = 0 ), axis = 1 ) 
    clw = np.flip( np.flip( np.mean( constants.extract_data( f, 'clwc' ), axis = 0 ), axis = 0 ), axis = 1 )  # kg/kg
    cli = np.flip( np.flip( np.mean( constants.extract_data( f, 'ciwc' ), axis = 0 ), axis = 0 ), axis = 1 )  # kg/kg
    ta = np.flip( np.flip( np.mean( constants.extract_data( f, 't' ), axis = 0), axis = 0 ), axis = 1 ) 
    w = np.flip( np.flip( np.mean( constants.extract_data( f, 'w' ), axis = 0), axis = 0 ), axis = 1 )  # Pa/s


cl_alt_lat = np.mean( cl, axis = -1 ) 
full_clw_alt_lat = np.mean( clw, axis = -1 )  # kg/kg
cli_alt_lat = np.mean( cli, axis = -1 )  # kg/kg
full_ta_alt_lat = np.mean( ta, axis = -1 ) 
w_alt_lat = np.mean( w, axis = -1 )  # Pa/s

#---create southern ocean and global datasets---#
    
cl_g = constants.global3DMean(cl, raw_lat)
cl_so = constants.global3DMean(cl[:,start_idx:end_idx], raw_lat[start_idx:end_idx])

clw_g = constants.global3DMean(clw, raw_lat)
clw_so = constants.global3DMean(clw[:,start_idx:end_idx], raw_lat[start_idx:end_idx])

cli_g = constants.global3DMean(cli, raw_lat)
cli_so = constants.global3DMean(cli[:,start_idx:end_idx], raw_lat[start_idx:end_idx])


#-------------calculate density---------------#

    ######## Binned Temperature Data ########
# values to bin: clw_alt_lat and ta_alt_lat
# binned into constants.ta_g array
# values in each bin to be summed
# call the summed values clw_t_g

stat = 'mean'
clw_t_g, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat.flatten(), full_clw_alt_lat.flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))
clw_t_so, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat[:,start_idx:end_idx].flatten(), full_clw_alt_lat[:,start_idx:end_idx].flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))

fig, ax = plt.subplots()
ax.plot( constants.ta, clw_t_g )
ax.plot( constants.ta, clw_t_so )
ax.axvline(x=273, label = '273K', color = 'black', linestyle='--')
plt.grid(True)
plt.show()
######################


#---interpolate and fit data---#

interpolated = interpolate.interp1d(raw_alt, cl_g, kind = 'cubic', fill_value="extrapolate")
cl_g = interpolated(constants.alt)

interpolated = interpolate.interp1d(raw_alt, clw_g, kind = 'cubic', fill_value="extrapolate")
clw_g = interpolated(constants.liq_alt)

interpolated = interpolate.interp1d(raw_alt, cli_g, kind = 'cubic', fill_value="extrapolate")
cli_g = interpolated(constants.alt)

interpolated = interpolate.interp1d(raw_alt, cl_so, kind = 'cubic', fill_value="extrapolate")
cl_so = interpolated(constants.alt)

interpolated = interpolate.interp1d(raw_alt, clw_so, kind = 'cubic', fill_value="extrapolate")
clw_so = interpolated(constants.liq_alt)

interpolated = interpolate.interp1d(raw_alt, cli_so, kind = 'cubic', fill_value="extrapolate")
cli_so = interpolated(constants.alt)


interpolated = interpolate.interp2d( raw_alt, raw_lat, np.transpose( cl_alt_lat ), kind = 'cubic')
cl_alt_lat = interpolated( constants.alt, constants.lat )
cl_alt_lat[ cl_alt_lat < 0 ] = np.nan

interpolated = interpolate.interp2d( raw_alt, raw_lat, np.transpose( full_clw_alt_lat ), kind = 'cubic')
clw_alt_lat = interpolated( constants.liq_alt, constants.lat )
clw_alt_lat[ clw_alt_lat < 0 ] = np.nan

interpolated = interpolate.interp2d( raw_alt, raw_lat, np.transpose( full_clw_alt_lat ), kind = 'cubic', fill_value="extrapolate")
full_clw_alt_lat = interpolated( constants.alt, constants.lat )
full_clw_alt_lat[ full_clw_alt_lat < 0 ] = np.nan

interpolated = interpolate.interp2d( raw_alt, raw_lat, np.transpose( cli_alt_lat ), kind = 'cubic')
cli_alt_lat = interpolated( constants.alt, constants.lat )
cli_alt_lat[ cli_alt_lat < 0 ] = np.nan

interpolated = interpolate.interp2d( raw_alt, raw_lat, np.transpose( full_ta_alt_lat ), kind = 'cubic', fill_value="extrapolate")
ta_alt_lat = interpolated( constants.liq_alt, constants.lat )
ta_alt_lat[ ta_alt_lat < 0 ] = np.nan

interpolated = interpolate.interp2d( raw_alt, raw_lat, np.transpose( full_ta_alt_lat ), kind = 'cubic', fill_value="extrapolate")
full_ta_alt_lat = interpolated( constants.alt, constants.lat )
full_ta_alt_lat[ full_ta_alt_lat < 0 ] = np.nan

interpolated = interpolate.interp2d( raw_alt, raw_lat, np.transpose( w_alt_lat ), kind = 'cubic', fill_value="extrapolate")
w_alt_lat = interpolated( constants.alt, constants.lat )

#----------------------------#



os.chdir( location + '/climate-analysis/reduced_data' )

save_filename = 'Jan_2006_Dec_2010_ECMWF.h5'

with h5py.File(save_filename, 'w') as p:

    p.create_dataset('clt', data=clt) # total cloud fraction corresponding to lat
    p.create_dataset('clt_lat_lon', data=clt_lat_lon ) # total cloud fraction corresponding to lat, lon
  
    p.create_dataset('clt_lc', data=clt_lc) # total cloud fraction corresponding to lat
    p.create_dataset('clt_lc_lat_lon', data=clt_lc_lat_lon ) # total cloud fraction corresponding to lat, lon

    p.create_dataset('clwvi', data=clwvi) # total cloud liquid water fraction corresponding to lat
    p.create_dataset('clwvi_lat_lon', data=clwvi_lat_lon ) # total cloud fraction corresponding to lat, lon

    p.create_dataset('clivi', data=clivi) # total cloud ice fraction corresponding to lat
    p.create_dataset('clivi_lat_lon', data=clivi_lat_lon ) # total cloud fraction corresponding to lat, lon
    
    p.create_dataset('cl_g', data=cl_g) # global layer total cloud fraction corresponding to alt
    p.create_dataset('clw_g', data=clw_g) # global layer cloud liquid water fraction corresponding to liq_alt
    p.create_dataset('cli_g', data=cli_g) # global layer cloud ice water fraction corresponding to alt
    
    p.create_dataset('cl_so', data=cl_so) # southern ocean layer total cloud fraction corresponding to alt
    p.create_dataset('clw_so', data=clw_so) # southern ocean layer cloud liquid water fraction corresponding to liq_alt
    p.create_dataset('cli_so', data=cli_so) # southern ocean layer cloud ice water fraction corresponding to alt

    p.create_dataset('clw_t_g', data=clw_t_g) # global layer cloud liquid water fraction corresponding to ta_g
    p.create_dataset('clw_t_so', data=clw_t_so) # global layer cloud liquid water fraction corresponding to ta_so
      
    p.create_dataset('w_alt_lat', data=np.transpose( w_alt_lat ) ) # vertical pressure velocity (Pa/s) corresponding to alt, lat

    p.create_dataset('ta_alt_lat', data= np.transpose( ta_alt_lat ) ) # temperature corresponding to liq_alt and lat
    p.create_dataset('cl_alt_lat', data=np.transpose( cl_alt_lat ) ) # total cloud fraction corresponding to alt and lat
    p.create_dataset('clw_alt_lat', data=np.transpose( clw_alt_lat ) ) # cloud liquid water fraction corresponding to liq_alt and lat
    p.create_dataset('cli_alt_lat', data=np.transpose( cli_alt_lat ) ) # cloud ice water fraction corresponding to alt and lat

    p.create_dataset('full_ta_alt_lat', data= np.transpose( full_ta_alt_lat ) ) # temperature corresponding to alt and lat
    p.create_dataset('full_clw_alt_lat', data=np.transpose( full_clw_alt_lat ) ) # total cloud fraction corresponding to alt and lat

    p.close()

end = time.time()
print('Create dataset took:', end - start, 's')

