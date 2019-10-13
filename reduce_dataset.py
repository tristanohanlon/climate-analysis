# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 18:40:37 2019

@author: Tristan O'Hanlon - University of Auckland & Jonathan Rogers

These functions will create a standard model dataset when parsed:
    
reduce_dataset( directory, filename, save_location, start, end, cl_is_fractional=False )

example:
reduce_dataset( 'gfdl_am4', '_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'E:/University/University/MSc/Models/climate-analysis/GFDL-AM4', datetime.datetime( 2000, 1, 1 ), datetime.datetime( 2006, 1, 1 ) )

if cl_is_fractional=True - the values get dived by 100
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


#---prime values for interpolation---#


def extract_data_over_time( type, f, start, end ):
        time_variable = f.variables['time']
        start_index = date2index( start, time_variable, select='before' )
        end_index = date2index( end, time_variable, select='before' )

        data = np.array( f.variables[type][start_index:end_index])
        return data

def extract_data( f, type ):
    return np.array( f.variables[type][:] )

def create_southern_ocean_data( raw_lat, global_data ):
    
    start_index = 0
    end_index = 0
    for index, l in enumerate( raw_lat ):
        if l < -70:
            start_index = index + 1
        if l <= -50:
            end_index = index
    so = global_data[start_index:end_index]
    return so


def reduce_dataset( directory, filename, save_location, start, end):
    os.chdir( directory )
  

    with Dataset( 'clt' + filename, 'r') as f:
        raw_lat = extract_data( f, 'lat')
        raw_lon = extract_data( f, 'lon')
        clt = extract_data_over_time('clt', f, start, end )
        if directory == 'CMIP6-GFDL-AM4':
            pass
        else:
            clt = clt / 100

    with Dataset( 'clwvi' + filename, 'r') as f:
        clwvi = extract_data_over_time('clwvi', f, start, end )

    with Dataset( 'clivi' + filename, 'r') as f:
        clivi = extract_data_over_time('clivi', f, start, end )


    with Dataset( 'cl' + filename, 'r') as f:
        cl = extract_data_over_time( 'cl',f, start, end )
        cl = np.mean( cl, axis = 0 ) # average over time
        if directory == 'CMIP5-CESM1-CAM5':
            pass
        else:
            cl = cl / 100
            
        if directory == 'CMIP6-GFDL-AM4' or directory == 'CMIP5-IPSL-CM5A-LR' or directory == 'CMIP6-IPSL-CM6A-LR':
            a = extract_data( f, 'ap')
            b = extract_data( f, 'b')
        else:
            a = extract_data( f, 'a')
            b = extract_data( f, 'b')

        if directory == 'CMIP6-GFDL-AM4' or directory == 'CMIP5-IPSL-CM5A-LR' or directory == 'CMIP6-IPSL-CM6A-LR':
            pass
        else:
            p0 = np.array(f.variables['p0'][:]) #in hPa
    

    with Dataset( 'clw' + filename, 'r') as f:
        clw = np.nanmean( extract_data_over_time( 'clw', f, start, end ), axis = 0 ) # average over time

    with Dataset( 'cli' + filename, 'r') as f:
        cli = np.nanmean( extract_data_over_time( 'cli', f, start, end ), axis = 0 ) # average over time

    with Dataset( 'ps' + filename, 'r') as f:
        ps = np.mean( extract_data_over_time( 'ps', f, start, end ), axis = 0 ) # average over time

    with Dataset( 'ta' + filename, 'r' ) as f:
        ta = extract_data_over_time( 'ta', f, start, end )
        ta[ta > 400] = np.nan
        ta = np.nanmean(ta, axis = 0) # average over time  
        plev_t = extract_data( f, 'plev') / 100

    #---get SO surface pressure---#
    ps_so = create_southern_ocean_data( raw_lat, np.mean(ps, axis = -1)) # average over longitude
    ps_so = np.mean( ps_so, axis=0)
    ps = np.mean(ps, axis = 0) # average over latitude
    ps = np.mean(ps, axis = 0) # average to single global average

    if directory == 'CMIP6-GFDL-AM4' or directory == 'CMIP5-IPSL-CM5A-LR' or directory == 'CMIP6-IPSL-CM6A-LR':
       p = a + b*ps
       p = np.array(p / 100) #hPa
       p_so = a + b*ps_so
       p_so = np.array(p_so / 100) #hPa 
    else:
        p = a*p0 + b*ps
        p = np.array(p / 100) #hPa
        p_so = a*p0 + b*ps_so
        p_so = np.array(p_so / 100) #hPa
    
    #---Approximate missing nan values in temperature data---#  

  
   
    #---convert pressure levels to altitude---#
    #https://www.mide.com/pages/air-pressure-at-altitude-calculator
    #https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
    raw_alt = np.empty((p.size,1),dtype=float)
    state = 0
    if directory == 'CMIP6-CESM2-CAM6':
        p = np.flip(p, axis = 0)

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



    #-------------Create grid sets---------------#

    clt_lat_lon = np.mean(clt, axis = 0) # average over time
    clwvi_lat_lon = np.mean(clwvi, axis = 0) # average over time
    clivi_lat_lon = np.mean(clivi, axis = 0) # average over time

    if directory == 'CMIP6-CESM2-CAM6':
        cl_alt_lat = np.flip( np.nanmean( cl , axis = -1 ), axis = 0 ) # average over longitude
        clw_alt_lat = np.flip( np.nanmean( clw , axis = -1 ), axis = 0 ) # average over longitude
        cli_alt_lat = np.flip( np.nanmean( cli , axis = -1 ), axis = 0 ) # average over longitude      
    else:
        cl_alt_lat = np.nanmean(cl , axis = -1) # average over longitude
        clw_alt_lat = np.nanmean(clw , axis = -1) # average over longitude
        cli_alt_lat = np.nanmean(cli , axis = -1) # average over longitude      

    full_ta_alt_lat = np.nanmean(ta, axis = -1)

    lwp_frac_lat_lon = (clwvi_lat_lon / (clwvi_lat_lon + clivi_lat_lon)) * clt_lat_lon
    iwp_frac_lat_lon = (clivi_lat_lon / (clwvi_lat_lon + clivi_lat_lon)) * clt_lat_lon
    
    lw_frac_alt_lat = (clw_alt_lat / (clw_alt_lat + cli_alt_lat)) * cl_alt_lat
    iw_frac_alt_lat = (cli_alt_lat / (clw_alt_lat + cli_alt_lat)) * cl_alt_lat


    interpolated = interpolate.interp1d(p, raw_alt, fill_value="extrapolate", kind = 'cubic')
    alt_temp = interpolated(plev_t)

    if directory == 'CMIP6-IPSL-CM6A-LR':
        raw_alt = raw_alt[:79]
        interpolated = interpolate.interp1d(p[:40], raw_alt[:40], fill_value="extrapolate", kind = 'cubic')
        alt_temp = interpolated(plev_t)
    else:
        interpolated = interpolate.interp1d(p, raw_alt, fill_value="extrapolate", kind = 'cubic')
        alt_temp = interpolated(plev_t)

      #-------------create cloud fraction datasets---------------#

    clt = np.nanmean(clt_lat_lon , axis = -1)

    cl_alt_lat[cl_alt_lat < 0] = 0

    cl_g = np.nanmean(cl_alt_lat , axis = -1) # average over latitude
    cl_so = create_southern_ocean_data( raw_lat, np.transpose(cl_alt_lat) )
    cl_so = np.nanmean(cl_so , axis = 0)

    #-------------create liquid and ice cloud fraction datasets---------------#

    clwvi = np.nanmean(lwp_frac_lat_lon , axis = -1)
    clivi = np.nanmean(iwp_frac_lat_lon , axis = -1)

    lw_frac_alt_lat[lw_frac_alt_lat < 0] = 0

    clw_g = np.nanmean(lw_frac_alt_lat , axis = -1)  # average over latitude
    clw_so = create_southern_ocean_data( raw_lat, np.transpose(lw_frac_alt_lat) )
    clw_so = np.nanmean(clw_so , axis = 0)
    lw_frac_alt_lat[np.isnan(lw_frac_alt_lat)] = 0
    clw_g[np.isnan(clw_g)] = 0
    clw_so[np.isnan(clw_so)] = 0

    iw_frac_alt_lat[iw_frac_alt_lat < 0] = 0

    cli_g = np.nanmean(iw_frac_alt_lat , axis = -1)  # average over latitude
    cli_so = create_southern_ocean_data( raw_lat, np.transpose(iw_frac_alt_lat) )
    cli_so = np.nanmean(cli_so , axis = 0)
    iw_frac_alt_lat[np.isnan(iw_frac_alt_lat)] = 0
    cli_g[np.isnan(cli_g)] = 0
    cli_so[np.isnan(cli_so)] = 0
    
    #--------------temp datasets--------------#
    
    ta_liq_g = np.nanmean(full_ta_alt_lat , axis = -1) # average over latitude
    ta_liq_so = create_southern_ocean_data( raw_lat, np.transpose(full_ta_alt_lat) )
    ta_liq_so = np.nanmean(ta_liq_so , axis = 0)
  
   

    #-------------interpolate to common lat, alt and liq_alt---------------#

    interpolated = interpolate.interp1d(raw_lat, clt, kind = 'cubic')
    clt = interpolated(constants.lat)

    interpolated = interpolate.interp1d(raw_lat, clwvi, kind = 'cubic')
    clwvi = interpolated(constants.lat)

    interpolated = interpolate.interp1d(raw_lat, clivi, kind = 'cubic')
    clivi = interpolated(constants.lat)

    interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clt_lat_lon ), kind = 'cubic')
    clt_lat_lon = interpolated(constants.lat, constants.lon)

    interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( lwp_frac_lat_lon ), kind = 'cubic')
    clwvi_lat_lon = interpolated(constants.lat, constants.lon)

    interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( iwp_frac_lat_lon ), kind = 'cubic')
    clivi_lat_lon = interpolated(constants.lat, constants.lon)

    pprint(clt_lat_lon)
    pprint(clwvi_lat_lon)


    if directory == 'CMIP6-IPSL-CM6A-LR':
        cl_alt_lat = constants.fit_2d_data(np.transpose(cl_alt_lat), raw_lat, raw_alt)
        cl_alt_lat = constants.fill(cl_alt_lat)
    else:
        interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( cl_alt_lat ), kind = 'cubic')
        cl_alt_lat = interpolated(constants.alt, constants.lat)

    if directory == 'CMIP6-IPSL-CM6A-LR':
        interpolated = interpolate.interp2d(raw_alt[:40], raw_lat, np.transpose( lw_frac_alt_lat[:40] ), kind = 'cubic')
        full_clw_alt_lat = interpolated(constants.alt, constants.lat)
#    elif directory == 'CMIP6-GFDL-AM4':
#        full_clw_alt_lat = constants.fit_2d_data(np.transpose(lw_frac_alt_lat), raw_lat, raw_alt)
#        full_clw_alt_lat = constants.fill(full_clw_alt_lat)
    else:
        interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( lw_frac_alt_lat ), kind = 'cubic')
        full_clw_alt_lat = interpolated(constants.alt, constants.lat)

    if directory == 'CMIP6-IPSL-CM6A-LR':
        interpolated = interpolate.interp2d(raw_alt[:40], raw_lat, np.transpose( lw_frac_alt_lat[:40] ), kind = 'cubic')
        clw_alt_lat = interpolated(constants.liq_alt, constants.lat)
#    elif directory == 'CMIP6-GFDL-AM4':
#        clw_alt_lat = constants.fit_2d_liq_data(np.transpose(lw_frac_alt_lat), raw_lat, raw_alt)
#        clw_alt_lat = constants.fill(clw_alt_lat)
    else:
        interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( lw_frac_alt_lat ), kind = 'cubic')
        clw_alt_lat = interpolated(constants.liq_alt, constants.lat)

    if directory == 'CMIP6-IPSL-CM6A-LR':
        interpolated = interpolate.interp2d(raw_alt[:40], raw_lat, np.transpose( iw_frac_alt_lat[:40] ), kind = 'cubic')
        cli_alt_lat = interpolated(constants.alt, constants.lat)
#    elif directory == 'CMIP6-GFDL-AM4':
#        cli_alt_lat = constants.fit_2d_data(np.transpose(iw_frac_alt_lat), raw_lat, raw_alt)
#        cli_alt_lat = constants.fill(cli_alt_lat)
    else:    
        interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( iw_frac_alt_lat ), kind = 'cubic')
        cli_alt_lat = interpolated(constants.alt, constants.lat)

        
    imp = SimpleImputer(missing_values=np.nan, strategy='mean')
    imp.fit(np.transpose(full_ta_alt_lat))  
    a = imp.transform(np.transpose(full_ta_alt_lat))
    ta_fixed = np.transpose(a)    
  
    if directory == 'CMIP6-IPSL-CM6A-LR':
        interpolated = interpolate.interp2d(alt_temp, raw_lat, np.transpose( ta_fixed ), kind = 'cubic')
        full_ta_alt_lat = interpolated(constants.alt, constants.lat)
    else:
        interpolated = interpolate.interp2d(alt_temp, raw_lat, np.transpose( ta_fixed ), kind = 'cubic')
        full_ta_alt_lat = interpolated(constants.alt, constants.lat)

    if directory == 'CMIP6-IPSL-CM6A-LR':
        interpolated = interpolate.interp2d(alt_temp, raw_lat, np.transpose( ta_fixed ), kind = 'cubic')
        ta_alt_lat = interpolated(constants.liq_alt, constants.lat)
    else:
        interpolated = interpolate.interp2d(alt_temp, raw_lat, np.transpose( ta_fixed ), kind = 'cubic')
        ta_alt_lat = interpolated(constants.liq_alt, constants.lat)


    
#    ta_alt_lat = constants.fit_2d_liq_data(np.transpose(np.nanmean(ta, axis = -1)), raw_lat, alt_temp)
#    ta_alt_lat = constants.fill(ta_alt_lat)
#    
#    full_ta_alt_lat = constants.fit_2d_data(np.transpose(np.nanmean(ta, axis = -1)), raw_lat, alt_temp)
#    full_ta_alt_lat = constants.fill(full_ta_alt_lat)


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

    interpolated = interpolate.interp1d(alt_temp, ta_liq_g, kind = 'cubic', fill_value="extrapolate")
    ta_liq_g = interpolated(constants.liq_alt)

    interpolated = interpolate.interp1d(alt_temp, ta_liq_so, kind = 'cubic', fill_value="extrapolate")
    ta_liq_so = interpolated(constants.liq_alt)

  
    #----------------------------#
    interpolated = interpolate.interp1d(ta_liq_g, clw_g, kind = 'cubic', fill_value="extrapolate")
    clw_t_g = interpolated(constants.ta_g)
    clw_t_g[clw_t_g < 0] = np.nan

    interpolated = interpolate.interp1d(ta_liq_so, clw_so, kind = 'cubic', fill_value="extrapolate")
    clw_t_so = interpolated(constants.ta_so)
    clw_t_so[clw_t_so < 0] = np.nan

    #----------------------------#


    os.chdir(save_location)

    save_filename = start.strftime( '%b_%Y') + '_' + end.strftime( '%b_%Y') + '_' + directory + '.h5'

    with h5py.File(save_filename, 'w') as p:
        
        p.create_dataset('clt', data=clt) # total cloud fraction corresponding to lat
        p.create_dataset('clt_lat_lon', data=np.transpose( clt_lat_lon )) # total cloud fraction corresponding to lat, lon
  
        p.create_dataset('clwvi', data=clwvi) # total cloud liquid water fraction corresponding to lat
        p.create_dataset('clwvi_lat_lon', data=np.transpose( clwvi_lat_lon )) # total cloud fraction corresponding to lat, lon

        p.create_dataset('clivi', data=clivi) # total cloud ice fraction corresponding to lat
        p.create_dataset('clivi_lat_lon', data=np.transpose( clivi_lat_lon )) # total cloud fraction corresponding to lat, lon
      
        p.create_dataset('cl_g', data=cl_g) # global layer total cloud fraction corresponding to alt
        p.create_dataset('clw_g', data=clw_g) # global layer cloud liquid water fraction corresponding to liq_alt
        p.create_dataset('cli_g', data=cli_g) # global layer cloud ice water fraction corresponding to alt
        
        p.create_dataset('cl_so', data=cl_so) # southern ocean layer total cloud fraction corresponding to alt
        p.create_dataset('clw_so', data=clw_so) # southern ocean layer cloud liquid water fraction corresponding to liq_alt
        p.create_dataset('cli_so', data=cli_so) # southern ocean layer cloud ice water fraction corresponding to alt
 
        p.create_dataset('clw_t_g', data=clw_t_g) # global layer cloud liquid water fraction corresponding to ta_g
        p.create_dataset('clw_t_so', data=clw_t_so) # global layer cloud liquid water fraction corresponding to ta_so
       
        p.create_dataset('ta_alt_lat', data=np.transpose( ta_alt_lat )) # temperature corresponding to liq_alt and lat
        p.create_dataset('cl_alt_lat', data=np.transpose( cl_alt_lat ) ) # total cloud fraction corresponding to alt and lat
        p.create_dataset('clw_alt_lat', data=np.transpose( clw_alt_lat ) ) # cloud liquid water fraction corresponding to liq_alt and lat
        p.create_dataset('cli_alt_lat', data=np.transpose( cli_alt_lat ) ) # cloud ice water fraction corresponding to alt and lat

        p.create_dataset('full_ta_alt_lat', data= np.transpose( full_ta_alt_lat )) # temperature corresponding to alt and lat
        p.create_dataset('full_clw_alt_lat', data=np.transpose( full_clw_alt_lat ) ) # total cloud fraction corresponding to alt and lat
  
        p.close()
