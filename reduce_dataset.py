# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 18:40:37 2019

@author: Tristan O'Hanlon

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

    

    if directory == 'CMIP6-CESM2-CAM6':
        with Dataset( 'clt' + filename, 'r') as f:
            clt = extract_data_over_time( 'CLDTOT', f, start, end ) # average over time
            raw_lat = extract_data( f, 'lat')
            alt = constants.alt[::-1]
            liq_alt = constants.liq_alt[::-1]
    else:
        with Dataset( 'clt' + filename, 'r') as f:
            raw_lat = extract_data( f, 'lat')
            clt = extract_data_over_time('clt', f, start, end )

    clt = np.mean(clt, axis = 0) # average over time
    clt = np.mean(clt, axis = -1) # average over longitude


    if directory == 'CMIP6-CESM2-CAM6':
        with Dataset( 'cl' + filename, 'r') as f:
            cl = extract_data_over_time( 'CLOUD', f, start, end ) # average over time
            cl = np.mean( cl, axis = 0 ) # average over time
            a = extract_data( f, 'hyam')
            b = extract_data( f, 'hybm')      
            p0 = np.array(f.variables['P0'][:]) #in hPa

    else:
        with Dataset( 'cl' + filename, 'r') as f:
            cl = extract_data_over_time( 'cl',f, start, end )
            cl = np.mean( cl, axis = 0 ) # average over time
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
    

    if directory == 'CMIP6-CESM2-CAM6':
        with Dataset( 'clw' + filename, 'r') as f:
            clw = np.nanmean( extract_data_over_time( 'CLDLIQ', f, start, end ), axis = 0 ) # average over time
    else:
        with Dataset( 'clw' + filename, 'r') as f:
            clw = np.nanmean( extract_data_over_time( 'clw', f, start, end ), axis = 0 ) # average over time

    if directory == 'CMIP6-CESM2-CAM6':
        with Dataset( 'cli' + filename, 'r') as f:
            cli = np.nanmean( extract_data_over_time( 'CLDICE', f, start, end ), axis = 0 ) # average over time
    else:
        with Dataset( 'cli' + filename, 'r') as f:
            cli = np.nanmean( extract_data_over_time( 'cli', f, start, end ), axis = 0 ) # average over time

    if directory == 'CMIP6-CESM2-CAM6':
        with Dataset( 'ps' + filename, 'r') as f:
            ps = np.nanmean( extract_data_over_time( 'PS', f, start, end ), axis = 0 ) # average over time
    else:
        with Dataset( 'ps' + filename, 'r') as f:
            ps = np.mean( extract_data_over_time( 'ps', f, start, end ), axis = 0 ) # average over time

    if directory == 'CMIP6-CESM2-CAM6':
        with Dataset( 'ta' + filename, 'r') as f:
            ta = extract_data_over_time( 'T', f, start, end )
            ta[ta > 400] = np.nan
            ta = np.nanmean(ta, axis = 0) # average over time  
    else:
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

    imp = SimpleImputer(missing_values=np.nan, strategy='mean')
    imp.fit(np.transpose(np.nanmean(ta, axis = -1)))  
    a = imp.transform(np.transpose(np.nanmean(ta, axis = -1)))
    ta_fixed = np.transpose(a)    
  
    if directory == 'CMIP6-CESM2-CAM6':
        ta_fixed = np.flip( ta_fixed, axis = 0 ) # average over longitude
   
    #---convert pressure levels to altitude---#
    #https://www.mide.com/pages/air-pressure-at-altitude-calculator
    #https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
    raw_alt = np.empty((p.size,1),dtype=float)
    state = 0

    if directory == 'CMIP6-CESM2-CAM6':
        p = np.flip(p)

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

    if directory == 'CMIP6-IPSL-CM6A-LR':
        raw_alt = raw_alt[:40]
        p = p[:40]
        cl_alt_lat = np.nanmean(cl , axis = -1)[:40] # average over longitude
        clw_alt_lat = np.nanmean(clw , axis = -1)[:40] # average over longitude
        cli_alt_lat = np.nanmean(cli , axis = -1)[:40] # average over longitude      
    if directory == 'CMIP6-CESM2-CAM6':
        cl_alt_lat = np.flip( np.nanmean( cl , axis = -1 ), axis = 0 ) # average over longitude
        clw_alt_lat = np.flip( np.nanmean( clw , axis = -1 ), axis = 0 ) # average over longitude
        cli_alt_lat = np.flip( np.nanmean( cli , axis = -1 ), axis = 0 ) # average over longitude      
    else:
        cl_alt_lat = np.nanmean(cl , axis = -1) # average over longitude
        clw_alt_lat = np.nanmean(clw , axis = -1) # average over longitude
        cli_alt_lat = np.nanmean(cli , axis = -1) # average over longitude      
    
    #interpolate corresponding temp altitudes

    if directory == 'CMIP6-CESM2-CAM6':
        alt_temp = raw_alt
    else:
        interpolated = interpolate.interp1d(p, raw_alt, fill_value="extrapolate", kind = 'cubic')
        alt_temp = interpolated(plev_t)
   

    #-------------interpolate to common lat, alt and liq_alt---------------#

    interpolated = interpolate.interp1d(raw_lat, clt, kind = 'cubic', fill_value="extrapolate")
    clt = interpolated(constants.lat)

    interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( cl_alt_lat ), kind = 'cubic')
    cl_alt_lat = interpolated(constants.alt, constants.lat)

    interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( clw_alt_lat ), kind = 'cubic')
    clw_alt_lat = interpolated(constants.alt, constants.lat)
    
    interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( cli_alt_lat ), kind = 'cubic')
    cli_alt_lat = interpolated(constants.alt, constants.lat)

    if directory == 'CMIP6-IPSL-CM6A-LR':
        interpolated = interpolate.interp2d(alt_temp, raw_lat, np.transpose( ta_fixed ), kind = 'linear')
        ta_alt_lat = interpolated(constants.liq_alt, constants.lat)
    else:
        interpolated = interpolate.interp2d(alt_temp, raw_lat, np.transpose( ta_fixed ), kind = 'cubic')
        ta_alt_lat = interpolated(constants.liq_alt, constants.lat)
        
    
    
    lw_frac_alt_lat = (clw_alt_lat / (clw_alt_lat + cli_alt_lat)) * cl_alt_lat
    iw_frac_alt_lat = (cli_alt_lat / (clw_alt_lat + cli_alt_lat)) * cl_alt_lat

    interpolated = interpolate.interp2d(constants.alt, constants.lat, lw_frac_alt_lat, kind = 'cubic')
    clw_alt_lat = interpolated(constants.liq_alt, constants.lat)


    #-------------create cloud fraction datasets---------------#
    cl_alt_lat[cl_alt_lat < 0] = 0

    cl_g = np.nanmean(cl_alt_lat , axis = 0) # average over latitude
    cl_so = create_southern_ocean_data( constants.lat, cl_alt_lat )
    cl_so = np.nanmean(cl_so , axis = 0)

    #-------------create liquid and ice cloud fraction datasets---------------#
    clw_alt_lat[clw_alt_lat < 0] = 0

    clw_g = np.nanmean(clw_alt_lat , axis = 0)  # average over latitude
    clw_so = create_southern_ocean_data( constants.lat, clw_alt_lat )
    clw_so = np.nanmean(clw_so , axis = 0)


    iw_frac_alt_lat[iw_frac_alt_lat < 0] = 0

    cli_g = np.nanmean(iw_frac_alt_lat , axis = 0)  # average over latitude
    cli_so = create_southern_ocean_data( constants.lat, iw_frac_alt_lat )
    cli_so = np.nanmean(cli_so , axis = 0)
     
    #--------------temp datasets--------------#
    
    ta_liq_g = np.nanmean(ta_alt_lat , axis = 0) # average over latitude
    ta_liq_so = create_southern_ocean_data( constants.lat, ta_alt_lat )
    ta_liq_so = np.nanmean(ta_liq_so , axis = 0)
   
    #----------------------------#

    os.chdir(save_location)

    save_filename = start.strftime( '%b_%Y') + '_' + end.strftime( '%b_%Y') + '_' + directory + '.h5'

    with h5py.File(save_filename, 'w') as p:
        
        p.create_dataset('ta_liq_g', data=ta_liq_g) # global layer temperature corresponding to liq_alt
        p.create_dataset('ta_liq_so', data=ta_liq_so) # southern ocean layer temperature corresponding to liq_alt
           
        p.create_dataset('clt', data=clt) # total cloud fraction corresponding to lat
        
        p.create_dataset('cl_g', data=cl_g) # global layer total cloud fraction corresponding to alt
        p.create_dataset('clw_g', data=clw_g) # global layer cloud liquid water fraction corresponding to liq_alt
        p.create_dataset('cli_g', data=cli_g) # global layer cloud ice water fraction corresponding to alt
        
        p.create_dataset('cl_so', data=cl_so) # southern ocean layer total cloud fraction corresponding to alt
        p.create_dataset('clw_so', data=clw_so) # southern ocean layer cloud liquid water fraction corresponding to liq_alt
        p.create_dataset('cli_so', data=cli_so) # southern ocean layer cloud ice water fraction corresponding to alt
        
        p.create_dataset('ta_alt_lat', data= np.transpose( ta_alt_lat )) # temperature corresponding to alt and lat
        p.create_dataset('cl_alt_lat', data=np.transpose( cl_alt_lat ) ) # total cloud fraction corresponding to alt and lat
        p.create_dataset('clw_alt_lat', data=np.transpose( clw_alt_lat ) ) # cloud liquid water fraction corresponding to liq_alt and lat
        p.create_dataset('cli_alt_lat', data=np.transpose( iw_frac_alt_lat ) ) # cloud ice water fraction corresponding to alt and lat
    
        p.close()
