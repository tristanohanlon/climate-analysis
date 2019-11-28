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
  
    ################################################---get regional data (time, lat, lon)---################################################


    with Dataset( 'clt' + filename, 'r') as f:
        raw_lat = extract_data( f, 'lat')
        raw_lon = extract_data( f, 'lon')
        clt_lat_lon = np.mean( extract_data_over_time('clt', f, start, end ), axis = 0 ) # average over time
        if directory == 'CMIP6-GFDL-AM4':
            pass
        else:
            clt_lat_lon = clt_lat_lon / 100

    with Dataset( 'clwvi' + filename, 'r') as f:
        clwvi = np.mean( extract_data_over_time('clwvi', f, start, end ), axis = 0 ) # average over time

    with Dataset( 'clivi' + filename, 'r') as f:
        clivi = np.mean( extract_data_over_time('clivi', f, start, end ), axis = 0 ) # average over time

    #---get radiative flux data---#

    with Dataset( 'rsdt' + filename, 'r') as f:
        rsdt = np.mean( extract_data_over_time('rsdt', f, start, end ), axis = 0 ) # average over time

    with Dataset( 'rsut' + filename, 'r') as f:
        rsut = np.mean( extract_data_over_time('rsut', f, start, end ), axis = 0 ) # average over time

    with Dataset( 'rsutcs' + filename, 'r') as f:
        rsutcs = np.mean( extract_data_over_time('rsutcs', f, start, end ), axis = 0 ) # average over time

    with Dataset( 'rlut' + filename, 'r') as f:
        rlut = np.mean( extract_data_over_time('rlut', f, start, end ), axis = 0 ) # average over time

    with Dataset( 'rlutcs' + filename, 'r') as f:
        rlutcs = np.mean( extract_data_over_time('rlutcs', f, start, end ), axis = 0 ) # average over time

    if directory == 'CMIP5-MRI-CGCM3' or directory == 'CMIP6-MRI-ESM2':
        rtmt = ( rsdt ) - ( rsut + rlut ) # https://www4.uwsp.edu/geo/faculty/lemke/geog101/lectures/02_radiation_energy_balance.html

    else:
        with Dataset( 'rtmt' + filename, 'r') as f:
            rtmt = np.mean( extract_data_over_time('rtmt', f, start, end ), axis = 0 ) # average over time


    #---create albedo data---#

    albedo_reg = rsut / rsdt
    # albedo_glob = np.nanmean( np.nanmean( albedo_reg, axis = 0 ), axis = 0 )
    albedo_so = create_southern_ocean_data( raw_lat, np.nanmean( albedo_reg, axis = 0 ) )
    # all_toa_net_glob = np.nanmean( np.nanmean( rtmt, axis = 0 ), axis = 0 )
    all_toa_net_so = create_southern_ocean_data( raw_lat, np.nanmean( rtmt, axis = 0 ) )
    # inc_solar = np.nanmean( np.nanmean( rsdt, axis = 0 ), axis = 0 )
    # ref_solar = np.nanmean( np.nanmean( rsut, axis = 0 ), axis = 0 )
    # out_LW = np.nanmean( np.nanmean( rlut, axis = 0 ), axis = 0 )

    # print( directory  )

    # print( str( inc_solar ) )
    # print( str( ref_solar ) )
    # print( str( out_LW ) )

    # print( str( all_toa_net_glob ) )
    # print( str( all_toa_net_so ) )
    # print( str( albedo_glob ) )
    # print( str( albedo_so ) )


    #---get aerosol data---#

    if directory == 'CMIP6-GFDL-AM4':
        with Dataset( 'mmrdust' + filename, 'r') as f:
            mmrdust = np.nanmean( extract_data_over_time( 'mmrdust', f, start, end ), axis = 0 ) # average over time
            mmrdust_lat_lon = np.nanmean(mmrdust, axis = 0) # average over altitude

        with Dataset( 'mmroa' + filename, 'r') as f:
            mmroa = np.nanmean( extract_data_over_time( 'mmroa', f, start, end ), axis = 0 ) # average over time
            mmroa_lat_lon = np.nanmean(mmroa, axis = 0) # average over altitude

        with Dataset( 'mmrso4' + filename, 'r') as f:
            mmrso4 = np.nanmean( extract_data_over_time( 'mmrso4', f, start, end ), axis = 0 ) # average over time
            mmrso4_lat_lon = np.nanmean(mmrso4, axis = 0) # average over altitude

        aerosol_norm = ( ( mmrdust_lat_lon / np.max( mmrdust_lat_lon ) ) + ( mmroa_lat_lon / np.max( mmroa_lat_lon ) ) + ( mmrso4_lat_lon / np.max( mmrso4_lat_lon ) ) ) / 3

    if directory == 'CMIP5-GISS-E2R' or directory == 'CMIP5-MIROC5':
        with Dataset( 'loaddust' + filename, 'r') as f:
            mmrdust_lat_lon = np.nanmean( extract_data_over_time( 'loaddust', f, start, end ), axis = 0 ) # average over time

        if directory == 'CMIP5-GISS-E2R':
            with Dataset( 'loadpoa' + filename, 'r') as f:
                mmroa_lat_lon = np.nanmean( extract_data_over_time( 'loadpoa', f, start, end ), axis = 0 ) # average over time
        else:
            with Dataset( 'loadoa' + filename, 'r') as f:
                mmroa_lat_lon = np.nanmean( extract_data_over_time( 'loadoa', f, start, end ), axis = 0 ) # average over time

        with Dataset( 'loadso4' + filename, 'r') as f:
            mmrso4_lat_lon = np.nanmean( extract_data_over_time( 'loadso4', f, start, end ), axis = 0 ) # average over time
 
        with Dataset( 'loadbc' + filename, 'r') as f:
            mmrbc_lat_lon = np.nanmean( extract_data_over_time( 'loadbc', f, start, end ), axis = 0 ) # average over time

        with Dataset( 'loadss' + filename, 'r') as f:
            mmrss_lat_lon = np.nanmean( extract_data_over_time( 'loadss', f, start, end ), axis = 0 ) # average over time

        aerosol_norm = ( ( mmrdust_lat_lon / np.max( mmrdust_lat_lon ) ) + ( mmroa_lat_lon / np.max( mmroa_lat_lon ) ) + ( mmrso4_lat_lon / np.max( mmrso4_lat_lon ) ) + ( mmrbc_lat_lon / np.max( mmrbc_lat_lon ) ) + ( mmrss_lat_lon / np.max( mmrss_lat_lon ) ) ) / 5


    #---reduce and interpolate regional data---#

    lwp_frac_lat_lon = (clwvi / (clwvi + clivi)) * clt_lat_lon
    iwp_frac_lat_lon = (clivi / (clwvi + clivi)) * clt_lat_lon

    clt = np.nanmean(clt_lat_lon , axis = -1)
    clwvi = np.nanmean(lwp_frac_lat_lon , axis = -1)
    clivi = np.nanmean(iwp_frac_lat_lon , axis = -1)

    interpolated = interpolate.interp1d(raw_lat, clt, kind = 'cubic', fill_value = 'extrapolate')
    clt = interpolated(constants.lat)

    interpolated = interpolate.interp1d(raw_lat, clwvi, kind = 'cubic', fill_value = 'extrapolate')
    clwvi = interpolated(constants.lat)

    interpolated = interpolate.interp1d(raw_lat, clivi, kind = 'cubic', fill_value = 'extrapolate')
    clivi = interpolated(constants.lat)

    if directory == 'CMIP5-IPSL-CM5A-LR' or directory == 'CMIP6-IPSL-CM6A-LR':
        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clt_lat_lon ), kind = 'cubic', fill_value = None)
        clt_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( lwp_frac_lat_lon ), kind = 'cubic', fill_value = None)
        clwvi_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( iwp_frac_lat_lon ), kind = 'cubic', fill_value = None)
        clivi_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( rtmt ), kind = 'cubic', fill_value = None)
        rtmt_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( rsdt ), kind = 'cubic', fill_value = None)
        rsdt_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( rsut ), kind = 'cubic', fill_value = None)
        rsut_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( rsutcs ), kind = 'cubic', fill_value = None)
        rsutcs_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( albedo_reg ), kind = 'cubic', fill_value = None)
        albedo_reg = interpolated(constants.lat, constants.lon)

    else:
        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clt_lat_lon ), kind = 'cubic', fill_value = np.nan)
        clt_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( lwp_frac_lat_lon ), kind = 'cubic', fill_value = np.nan)
        clwvi_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( iwp_frac_lat_lon ), kind = 'cubic', fill_value = np.nan)
        clivi_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( rtmt ), kind = 'cubic', fill_value = np.nan)
        rtmt_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( rsdt ), kind = 'cubic', fill_value = np.nan)
        rsdt_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( rsut ), kind = 'cubic', fill_value = np.nan)
        rsut_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( rsutcs ), kind = 'cubic', fill_value = np.nan)
        rsutcs_lat_lon = interpolated(constants.lat, constants.lon)
        
        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( albedo_reg ), kind = 'cubic', fill_value = np.nan)
        albedo_reg = interpolated(constants.lat, constants.lon)

    if directory == 'CMIP5-MRI-CGCM3' or directory == 'CMIP6-MRI-ESM2':
        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( rlut ), kind = 'cubic', fill_value = np.nan)
        rlut_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( rlutcs ), kind = 'cubic', fill_value = np.nan)
        rlutcs_lat_lon = interpolated(constants.lat, constants.lon)


    if directory == 'CMIP6-GFDL-AM4':
        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( mmrdust_lat_lon ), kind = 'cubic', fill_value = np.nan)
        mmrdust_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( mmroa_lat_lon ), kind = 'cubic', fill_value = np.nan)
        mmroa_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( mmrso4_lat_lon ), kind = 'cubic', fill_value = np.nan)
        mmrso4_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( aerosol_norm ), kind = 'cubic', fill_value = np.nan)
        aerosol_norm_lat_lon = interpolated(constants.lat, constants.lon)

    if directory == 'CMIP5-GISS-E2R' or directory == 'CMIP5-MIROC5':
        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( mmrdust_lat_lon ), kind = 'cubic', fill_value = np.nan)
        mmrdust_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( mmroa_lat_lon ), kind = 'cubic', fill_value = np.nan)
        mmroa_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( mmrso4_lat_lon ), kind = 'cubic', fill_value = np.nan)
        mmrso4_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( mmrbc_lat_lon ), kind = 'cubic', fill_value = np.nan)
        mmrbc_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( mmrss_lat_lon ), kind = 'cubic', fill_value = np.nan)
        mmrss_lat_lon = interpolated(constants.lat, constants.lon)

        interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( aerosol_norm ), kind = 'cubic', fill_value = np.nan)
        aerosol_norm_lat_lon = interpolated(constants.lat, constants.lon)




    ################################################---get profile data (time, plevel, lat, lon)---################################################

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


    #---get global and so surface pressure---#

    ps_so = create_southern_ocean_data( raw_lat, np.mean(ps, axis = -1)) # average over longitude
    ps_so = np.mean( ps_so, axis=0)
    ps = np.mean(ps, axis = 0) # average over latitude
    ps = np.mean(ps, axis = 0) # average to single global average


    #---convert hybrid pressure values into hPa---#

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

    # get index for pressure just below 700hPa to define low cloud boundary
    for p_index, a in enumerate (p):
        if a >= 700:
            index = p_index
    

    raw_alt = np.transpose( raw_alt )[0]

    #-------------create alt_lat grid sets---------------#

    if directory == 'CMIP6-CESM2-CAM6':
        cl_alt_lat = np.flip( np.nanmean( cl , axis = -1 ), axis = 0 ) # average over longitude
        clw_alt_lat = np.flip( np.nanmean( clw , axis = -1 ), axis = 0 ) # average over longitude
        cli_alt_lat = np.flip( np.nanmean( cli , axis = -1 ), axis = 0 ) # average over longitude    
        cl = np.flip( cl, axis = 0 )  
        clw = np.flip( clw, axis = 0 )  

    else:
        cl_alt_lat = np.nanmean(cl , axis = -1) # average over longitude
        clw_alt_lat = np.nanmean(clw , axis = -1) # average over longitude
        cli_alt_lat = np.nanmean(cli , axis = -1) # average over longitude      

    full_ta_alt_lat = np.nanmean(ta, axis = -1)

    
    #-------------create temp_alt grid sets---------------#


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

    # define low cloud cover (below 3km)
    clt_lc = cl[0:index]
    clt_lc_lat_lon = np.nanmean( clt_lc, axis = 0 )
    clt_lc = np.nanmean( clt_lc_lat_lon, axis = -1 )

    cl_alt_lat[cl_alt_lat < 0] = 0

    cl_g = np.nanmean(cl_alt_lat , axis = -1) # average over latitude
    cl_so = create_southern_ocean_data( raw_lat, np.transpose(cl_alt_lat) )
    cl_so = np.nanmean(cl_so , axis = 0)




    #-------------create liquid and ice cloud fraction datasets---------------#

    clw_lc = clw[0:index]
    clw_lc_lat_lon = np.nanmean( clw_lc, axis = 0 )

    cli_lc = cli[0:index]
    cli_lc_lat_lon = np.nanmean( cli_lc, axis = 0 )

    clwvi_lc_lat_lon = (clw_lc_lat_lon / (clw_lc_lat_lon + cli_lc_lat_lon)) * clt_lc_lat_lon
    clwvi_lc = np.nanmean( clwvi_lc_lat_lon, axis = -1 )

    if directory == 'CMIP6-MIROC6':
        clwvi_lc_lat_lon = constants.fill(clwvi_lc_lat_lon)

    lw_frac_alt_lat = (clw_alt_lat / (clw_alt_lat + cli_alt_lat)) * cl_alt_lat
    iw_frac_alt_lat = (cli_alt_lat / (clw_alt_lat + cli_alt_lat)) * cl_alt_lat

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
    
    
    #-------------interpolate liquid and ice fractions to common lat, alt and liq_alt---------------#


    interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clt_lc_lat_lon ), kind = 'cubic')
    clt_lc_lat_lon = interpolated(constants.lat, constants.lon)

    interpolated = interpolate.interp2d(raw_lat, raw_lon, np.transpose( clwvi_lc_lat_lon ), kind = 'cubic')
    clwvi_lc_lat_lon = interpolated(constants.lat, constants.lon)

    if directory == 'CMIP6-IPSL-CM6A-LR':
        cl_alt_lat = constants.fit_2d_data(np.transpose(cl_alt_lat), raw_lat, raw_alt)
        cl_alt_lat = constants.fill(cl_alt_lat)
    else:
        interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( cl_alt_lat ), kind = 'cubic')
        cl_alt_lat = interpolated(constants.alt, constants.lat)

    if directory == 'CMIP6-IPSL-CM6A-LR':
        interpolated = interpolate.interp2d(raw_alt[:40], raw_lat, np.transpose( lw_frac_alt_lat[:40] ), kind = 'cubic')
        full_clw_alt_lat = interpolated(constants.alt, constants.lat)
    else:
        interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( lw_frac_alt_lat ), kind = 'cubic')
        full_clw_alt_lat = interpolated(constants.alt, constants.lat)

    if directory == 'CMIP6-IPSL-CM6A-LR':
        interpolated = interpolate.interp2d(raw_alt[:40], raw_lat, np.transpose( lw_frac_alt_lat[:40] ), kind = 'cubic')
        clw_alt_lat = interpolated(constants.liq_alt, constants.lat)
    else:
        interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( lw_frac_alt_lat ), kind = 'cubic')
        clw_alt_lat = interpolated(constants.liq_alt, constants.lat)

    if directory == 'CMIP6-IPSL-CM6A-LR':
        interpolated = interpolate.interp2d(raw_alt[:40], raw_lat, np.transpose( iw_frac_alt_lat[:40] ), kind = 'cubic')
        cli_alt_lat = interpolated(constants.alt, constants.lat)
    else:    
        interpolated = interpolate.interp2d(raw_alt, raw_lat, np.transpose( iw_frac_alt_lat ), kind = 'cubic')
        cli_alt_lat = interpolated(constants.alt, constants.lat)

    interpolated = interpolate.interp1d(raw_lat, clt_lc, kind = 'cubic', fill_value="extrapolate")
    clt_lc = interpolated(constants.lat)

    interpolated = interpolate.interp1d(raw_lat, clwvi_lc, kind = 'cubic', fill_value="extrapolate")
    clwvi_lc = interpolated(constants.lat)

    interpolated = interpolate.interp1d(raw_alt, cl_g, kind = 'cubic', fill_value="extrapolate")
    cl_g = interpolated(constants.alt)
    if directory == 'CMIP6-IPSL-CM6A-LR' or directory == 'CMIP5-IPSL-CM5A-LR' or directory == 'CMIP6-MRI-ESM2':
        cl_g[:2] = np.nan
 

    interpolated = interpolate.interp1d(raw_alt, clw_g, kind = 'cubic', fill_value="extrapolate")
    clw_g = interpolated(constants.liq_alt)
    clw_g[clw_g < 0] = np.nan

    interpolated = interpolate.interp1d(raw_alt, cli_g, kind = 'cubic', fill_value="extrapolate")
    cli_g = interpolated(constants.alt)
    cli_g[cli_g < 0] = np.nan
    cli_g[:2] = np.nan
    
    interpolated = interpolate.interp1d(raw_alt, cl_so, kind = 'cubic', fill_value="extrapolate")
    cl_so = interpolated(constants.alt)
    if directory == 'CMIP6-IPSL-CM6A-LR' or directory == 'CMIP5-IPSL-CM5A-LR':
        cl_so[:2] = np.nan

    interpolated = interpolate.interp1d(raw_alt, clw_so, kind = 'cubic', fill_value="extrapolate")
    clw_so = interpolated(constants.liq_alt)
    clw_so[clw_so < 0] = np.nan

    interpolated = interpolate.interp1d(raw_alt, cli_so, kind = 'cubic', fill_value="extrapolate")
    cli_so = interpolated(constants.alt)
    cli_so[cli_so < 0] = np.nan
    cli_so[:2] = np.nan


    #--------------temp datasets--------------#
    
    ta_liq_g = np.nanmean(full_ta_alt_lat , axis = -1) # average over latitude
    ta_liq_so = create_southern_ocean_data( raw_lat, np.transpose(full_ta_alt_lat) )
    ta_liq_so = np.nanmean(ta_liq_so , axis = 0) # average over latitude
  

    #---Approximate missing nan values in temperature data---#  
       
    imp = SimpleImputer(missing_values=np.nan, strategy='mean')
    imp.fit(np.transpose(full_ta_alt_lat))  
    a = imp.transform(np.transpose(full_ta_alt_lat))
    ta_fixed = np.transpose(a)    
  

    #-------------interpolate temp data to common lat, alt and liq_alt---------------#

    interpolated = interpolate.interp2d(alt_temp, raw_lat, np.transpose( ta_fixed ), kind = 'cubic')
    full_ta_alt_lat = interpolated(constants.alt, constants.lat)

    interpolated = interpolate.interp2d(alt_temp, raw_lat, np.transpose( ta_fixed ), kind = 'cubic')
    ta_alt_lat = interpolated(constants.liq_alt, constants.lat)

    interpolated = interpolate.interp1d(alt_temp, ta_liq_g, kind = 'cubic', fill_value="extrapolate")
    ta_liq_g = interpolated(constants.liq_alt)

    interpolated = interpolate.interp1d(alt_temp, ta_liq_so, kind = 'cubic', fill_value="extrapolate")
    ta_liq_so = interpolated(constants.liq_alt)

  
    #----------fit clw to temp levels------------#

    interpolated = interpolate.interp1d(ta_liq_g, clw_g, kind = 'cubic', fill_value="extrapolate")
    clw_t_g = interpolated(constants.ta_g)
    clw_t_g[clw_t_g < 0] = np.nan

    clw_t_so = constants.fit_ta_so_data( clw_so, ta_liq_so )  
    clw_t_so[clw_t_so < 0] = np.nan
    clw_t_so_filled = constants.fill( clw_t_so[:36] )
    clw_t_so[:36] = clw_t_so_filled
    
    #----------------------------#


    os.chdir(save_location)

    save_filename = start.strftime( '%b_%Y') + '_' + end.strftime( '%b_%Y') + '_' + directory + '.h5'

    with h5py.File(save_filename, 'w') as p:
        
        p.create_dataset('clt', data=clt) # total cloud fraction corresponding to lat
        p.create_dataset('clt_lc', data=clt_lc) # low cloud fraction corresponding to lat
        p.create_dataset('clt_lat_lon', data=np.transpose( clt_lat_lon )) # total cloud fraction corresponding to lat, lon
        p.create_dataset('clt_lc_lat_lon', data=np.transpose( clt_lc_lat_lon )) # low cloud fraction corresponding to lat, lon
 
        p.create_dataset('clwvi', data=clwvi) # total cloud liquid water fraction corresponding to lat
        p.create_dataset('clwvi_lc', data=np.transpose( clwvi_lc )) # low cloud liquid water fraction corresponding to lat
        p.create_dataset('clwvi_lat_lon', data=np.transpose( clwvi_lat_lon )) # total cloud liquid water fraction corresponding to lat, lon
        p.create_dataset('clwvi_lc_lat_lon', data=np.transpose( clwvi_lc_lat_lon )) # low cloud liquid water fraction corresponding to lat, lon

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
        
        p.create_dataset('rsdt_lat_lon', data=np.transpose( rsdt_lat_lon )) # toa incoming sw flux corresponding to lat, lon
        p.create_dataset('rtmt_lat_lon', data=np.transpose( rtmt_lat_lon ) ) # net downwards raditaive flux at toa corresponding to lat, lon
        p.create_dataset('rsut_lat_lon', data=np.transpose( rsut_lat_lon ) ) # toa outgoing sw flux corresponding to lat, lon
        p.create_dataset('rsutcs_lat_lon', data=np.transpose( rsutcs_lat_lon ) ) # toa outgoing sw flux assuming clear sky corresponding to lat, lon
        p.create_dataset('albedo_reg', data=np.transpose( albedo_reg ) ) 
        p.create_dataset('all_toa_net_so', data=all_toa_net_so) 
        p.create_dataset('albedo_so', data=albedo_so) 

        if directory == 'CMIP5-MRI-CGCM3' or directory == 'CMIP6-MRI-ESM2':
            p.create_dataset('rlut_lat_lon', data=np.transpose( rlut_lat_lon ) ) # toa outgoing sw flux corresponding to lat, lon
            p.create_dataset('rlutcs_lat_lon', data=np.transpose( rlutcs_lat_lon ) ) # toa outgoing sw flux assuming clear sky corresponding to lat, lon

        if directory == 'CMIP6-GFDL-AM4':
            p.create_dataset('mmrdust_lat_lon', data=np.transpose( mmrdust_lat_lon )) # aerosol dust mixing ratio corresponding to lat, lon
            p.create_dataset('mmroa_lat_lon', data=np.transpose( mmroa_lat_lon ) ) # aerosol total organic mixing ratio corresponding to lat, lon
            p.create_dataset('mmrso4_lat_lon', data=np.transpose( mmrso4_lat_lon ) ) # aerosol so4 mixing ratio corresponding to lat, lon
            p.create_dataset('aerosol_norm_lat_lon', data=np.transpose( aerosol_norm_lat_lon ) ) # total normalised aerosols corresponding to lat, lon

        if directory == 'CMIP5-GISS-E2R' or directory == 'CMIP5-MIROC5':
            p.create_dataset('mmrdust_lat_lon', data=np.transpose( mmrdust_lat_lon )) # aerosol dust mixing ratio corresponding to lat, lon
            p.create_dataset('mmroa_lat_lon', data=np.transpose( mmroa_lat_lon ) ) # aerosol total organic mixing ratio corresponding to lat, lon
            p.create_dataset('mmrso4_lat_lon', data=np.transpose( mmrso4_lat_lon ) ) # aerosol so4 mixing ratio corresponding to lat, lon
            p.create_dataset('mmrbc_lat_lon', data=np.transpose( mmrbc_lat_lon ) ) # aerosol total organic mixing ratio corresponding to lat, lon
            p.create_dataset('mmrss_lat_lon', data=np.transpose( mmrss_lat_lon ) ) # aerosol so4 mixing ratio corresponding to lat, lon
            p.create_dataset('aerosol_norm_lat_lon', data=np.transpose( aerosol_norm_lat_lon ) ) # total normalised aerosols corresponding to lat, lon


        p.close()
