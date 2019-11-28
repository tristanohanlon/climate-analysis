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

lwp_frac_lat_lon = (clwvi_lat_lon / (clwvi_lat_lon + clivi_lat_lon)) * clt_lat_lon
iwp_frac_lat_lon = (clivi_lat_lon / (clwvi_lat_lon + clivi_lat_lon)) * clt_lat_lon
clwvi = np.nanmean(lwp_frac_lat_lon , axis = -1)
clivi = np.nanmean(iwp_frac_lat_lon , axis = -1)
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

interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( lwp_frac_lat_lon ), kind = 'cubic')
clwvi_lat_lon = interpolated(constants.lat, constants.lon)
clwvi_lat_lon = np.transpose(clwvi_lat_lon)

interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( iwp_frac_lat_lon ), kind = 'cubic')
clivi_lat_lon = interpolated(constants.lat, constants.lon)
clivi_lat_lon = np.transpose(clivi_lat_lon)


os.chdir( location + 'Data/ECMWF/pressure_levels' )
  

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


with Dataset( '200601-201012_ECMWF_plevel_T_cc_clw_cli_w.nc', 'r') as f:
    p = np.flip( extract_data( 'level', f ), axis = 0 )        
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
    raw_lat = np.flip( extract_data( 'latitude', f ), axis = 0 )


    cl_alt_lat = np.flip( np.flip( np.mean( np.mean( extract_data( 'cc', f ), axis = 0 ), axis = -1 ), axis = 0 ), axis = -1 ) 
    clw_alt_lat = np.flip( np.flip( np.mean( np.mean(extract_data( 'clwc', f ), axis = 0), axis = -1 ), axis = 0 ), axis = -1 )
    cli_alt_lat = np.flip( np.flip( np.mean( np.mean(extract_data( 'ciwc', f ), axis = 0), axis = -1 ), axis = 0 ), axis = -1 )
    full_ta_alt_lat = np.flip( np.flip( np.mean( np.mean(extract_data( 't', f ), axis = 0), axis = -1 ), axis = 0 ), axis = -1 )
    w_alt_lat = np.flip( np.flip( np.mean( np.mean(extract_data( 'w', f ), axis = 0), axis = -1 ), axis = 0 ), axis = -1 )

#---create liquid and ice fractions---#

full_clw_alt_lat = (clw_alt_lat / (clw_alt_lat + cli_alt_lat)) * cl_alt_lat
cli_alt_lat = (cli_alt_lat / (clw_alt_lat + cli_alt_lat)) * cl_alt_lat

#---interpolate and fit data---#

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

#---create southern ocean and global datasets---#
    
cl_so = create_southern_ocean_data( constants.lat, cl_alt_lat)
clw_so = create_southern_ocean_data( constants.lat, clw_alt_lat)
cli_so = create_southern_ocean_data( constants.lat, cli_alt_lat)
ta_liq_so = create_southern_ocean_data( constants.lat, ta_alt_lat)
w_so = create_southern_ocean_data( constants.lat, w_alt_lat)

cl_g = np.nanmean( cl_alt_lat, axis = 0 )
clw_g = np.nanmean( clw_alt_lat, axis = 0 )
cli_g = np.nanmean( cli_alt_lat, axis = 0 )
ta_liq_g = np.nanmean( ta_alt_lat, axis = 0 )
w_g = np.nanmean( w_alt_lat, axis = 0 )


#--------------create temp clw--------------#

interpolated = interpolate.interp1d(ta_liq_g, clw_g, kind = 'linear', fill_value="extrapolate")
clw_t_g = interpolated(constants.ta_g)
clw_t_g[clw_t_g < 0] = np.nan

interpolated = interpolate.interp1d(ta_liq_so, clw_so, kind = 'linear', fill_value="extrapolate")
clw_t_so = interpolated(constants.ta_so)
clw_t_so[clw_t_so < 0] = np.nan

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
    p.create_dataset('w_g', data=w_g) # global vertical pressure velocity (Pa/s) corresponding to alt
    p.create_dataset('w_so', data=w_so) # southern ocean vertical pressure velocity (Pa/s) corresponding to alt

    p.create_dataset('ta_alt_lat', data= np.transpose( ta_alt_lat ) ) # temperature corresponding to liq_alt and lat
    p.create_dataset('cl_alt_lat', data=np.transpose( cl_alt_lat ) ) # total cloud fraction corresponding to alt and lat
    p.create_dataset('clw_alt_lat', data=np.transpose( clw_alt_lat ) ) # cloud liquid water fraction corresponding to liq_alt and lat
    p.create_dataset('cli_alt_lat', data=np.transpose( cli_alt_lat ) ) # cloud ice water fraction corresponding to alt and lat

    p.create_dataset('full_ta_alt_lat', data= np.transpose( full_ta_alt_lat ) ) # temperature corresponding to alt and lat
    p.create_dataset('full_clw_alt_lat', data=np.transpose( full_clw_alt_lat ) ) # total cloud fraction corresponding to alt and lat

    p.close()

end = time.time()
print('Create dataset took:', end - start, 's')

