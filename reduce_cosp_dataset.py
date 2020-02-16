# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 18:40:37 2019

@author: Tristan O'Hanlon - University of Auckland & Jonathan Rogers

These functions will create a standard model dataset when parsed:
    
reduce_dataset( directory, filename, save_location, start, end, cl_is_fractional=False )

example:
reduce_dataset( 'gfdl_am4', '_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'E:/University/University/MSc/Models/climate-analysis/GFDL-AM4', datetime.datetime( 2000, 1, 1 ), datetime.datetime( 2006, 1, 1 ) )

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
from scipy import ndimage as nd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs



def reduce_cosp_dataset( directory, save_location, start, end, location ):

    os.chdir( location + 'Data/' + directory + '/COSP/' )
  
    ################################################---get regional data (time, lat, lon)---################################################

    print(directory)
    with Dataset( constants.variable_to_filename( 'cltcalipso' ), 'r') as f:
        raw_lat = constants.extract_data( 'lat', f )
        raw_lon = constants.extract_data( 'lon', f )

        # get indexes for southern ocean -70 to -50 lat

        start_idx = np.abs(raw_lat - (-69.5)).argmin()
        end_idx = np.abs(raw_lat - (-50)).argmin()

        raw_clt_lat_lon = np.mean( constants.extract_data_over_time('cltcalipso', f, start, end ), axis = 0 ) / 100 # average over time


    with Dataset( constants.variable_to_filename( 'clhcalipso' ), 'r') as f:
        raw_clt_h_lat_lon = np.mean( constants.extract_data_over_time('clhcalipso', f, start, end ), axis = 0 ) / 100 # average over time
        raw_clt_h_lat_lon[ raw_clt_h_lat_lon > 1 ] = np.nan

    with Dataset( constants.variable_to_filename( 'clmcalipso' ), 'r') as f:
        raw_clt_m_lat_lon = np.mean( constants.extract_data_over_time('clmcalipso', f, start, end ), axis = 0 ) / 100 # average over time
        raw_clt_m_lat_lon[ raw_clt_m_lat_lon > 1 ] = np.nan

    with Dataset( constants.variable_to_filename( 'cllcalipso' ), 'r') as f:
        raw_clt_l_lat_lon = np.mean( constants.extract_data_over_time('cllcalipso', f, start, end ), axis = 0 ) / 100 # average over time
        raw_clt_l_lat_lon[ raw_clt_l_lat_lon > 1 ] = np.nan

    imp = SimpleImputer(missing_values=np.nan, strategy='mean')
    imp.fit(np.transpose(raw_clt_l_lat_lon))  
    a = imp.transform(np.transpose(raw_clt_l_lat_lon))
    raw_clt_l_lat_lon = np.transpose(a)


    clt = np.mean( raw_clt_lat_lon, axis = -1 )
    clt_h = np.mean( raw_clt_h_lat_lon, axis = -1 )
    clt_m = np.mean( raw_clt_m_lat_lon, axis = -1 )
    clt_l = np.mean( raw_clt_l_lat_lon, axis = -1 )

    if raw_lat.shape[0] == constants.lat.shape[0]:
        pass
    else:
        interpolated = interpolate.interp1d(raw_lat, clt, kind = 'cubic', fill_value="extrapolate")
        clt = interpolated(constants.lat)
        interpolated = interpolate.interp1d(raw_lat, clt_h, kind = 'cubic', fill_value="extrapolate")
        clt_h = interpolated(constants.lat)
        interpolated = interpolate.interp1d(raw_lat, clt_m, kind = 'cubic', fill_value="extrapolate")
        clt_m = interpolated(constants.lat)
        if 'CAM' in directory or 'GISS-E2R' in directory:
            interpolated = interpolate.interp1d(raw_lat[1:], clt_l, kind = 'cubic', fill_value="extrapolate")
            clt_l = interpolated(constants.lat)
        else:
            interpolated = interpolate.interp1d(raw_lat, clt_l, kind = 'cubic', fill_value="extrapolate")
            clt_l = interpolated(constants.lat)

    interpolated = interpolate.interp2d(raw_lon, raw_lat, raw_clt_lat_lon, kind = 'cubic')
    clt_lat_lon = interpolated(constants.lon, constants.lat)

    interpolated = interpolate.interp2d(raw_lon, raw_lat, raw_clt_h_lat_lon, kind = 'cubic')
    clt_h_lat_lon = interpolated(constants.lon, constants.lat)

    interpolated = interpolate.interp2d(raw_lon, raw_lat, raw_clt_m_lat_lon, kind = 'cubic')
    clt_m_lat_lon = interpolated(constants.lon, constants.lat)
        
    if 'CAM' in directory or 'GISS-E2R' in directory:
        interpolated = interpolate.interp2d(raw_lon, raw_lat[1:], raw_clt_l_lat_lon, kind = 'linear')
        clt_l_lat_lon = interpolated(constants.lon, constants.lat)
    else:
        interpolated = interpolate.interp2d(raw_lon, raw_lat, raw_clt_l_lat_lon, kind = 'linear')
        clt_l_lat_lon = interpolated(constants.lon, constants.lat)


    # #----Test Plots----#


    # fig, ax = plt.subplots()
    # ax.plot( constants.lat[constants.lat_confine_1:constants.lat_confine_2], clt[constants.lat_confine_1:constants.lat_confine_2] )
    # ax.plot( constants.lat[constants.lat_confine_1:constants.lat_confine_2], clt_l[constants.lat_confine_1:constants.lat_confine_2] )

    # ax.set_ylabel('Cloud Fraction')
    # ax.set_xlabel('Latitude')
    # ax.set_title ('Global Cloud Fraction vs Latitude')
    # plt.grid(True)
    # plt.show()


    # ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    # ax.coastlines()
    # p = ax.contourf(constants.lon, constants.lat, clt_lat_lon, transform=ccrs.PlateCarree(), cmap='coolwarm')
    # cbar = plt.colorbar(p, orientation='horizontal')
    # cbar.set_label('Cloud Fraction')
    # ax.set_title('Total Cloud Fraction')
    # plt.show()  



    # import temperature data
    os.chdir( location + 'climate-analysis/reduced_data' )
    ta_data = start.strftime( '%b_%Y') + '_' + end.strftime( '%b_%Y') + '_' + directory + '.h5'
    h5f = h5py.File( ta_data, 'r')
    ta_alt_lat = h5f['ta_alt_lat'][:]
    full_ta_alt_lat = h5f['full_ta_alt_lat'][:]

    os.chdir( location + 'Data/' + directory + '/COSP/' )

    ######################---get profile data (time, alt40, lat, lon)---###################

    with Dataset( constants.variable_to_filename( 'clcalipso' ), 'r') as f:
        cl = constants.extract_data_over_time( 'clcalipso', f, start, end )
        cl = np.mean( cl, axis = 0 )  / 100 # average over time
        cl[ cl > 1 ] = np.nan

        if 'CM6A' in directory:
            raw_alt = constants.extract_data( 'height', f ) / 1000 # in km
        else:
            raw_alt = constants.extract_data( 'alt40', f ) / 1000 # in km


    cl_alt_lat = np.nanmean( cl, axis = -1 ) 

    imp = SimpleImputer(missing_values=np.nan, strategy='mean')
    imp.fit(np.transpose(cl_alt_lat))  
    a = imp.transform(np.transpose(cl_alt_lat))
    cl_alt_lat = np.transpose(a)    

    interpolated = interpolate.interp2d( raw_lat, raw_alt, cl_alt_lat, kind = 'cubic')
    cl_alt_lat = interpolated( constants.lat, raw_alt )
    cl_alt_lat[ cl_alt_lat <= 0 ] = 0

    cl_g = constants.globalalt_latMean(cl_alt_lat, constants.lat)
    cl_so = constants.globalalt_latMean(cl_alt_lat[:,start_idx:end_idx], constants.lat[start_idx:end_idx])


    # #----Test Plot----#

    # fig, ax = plt.subplots()
    # ax.plot( cl_so, raw_alt )
    # ax.plot( cl_g, raw_alt )
    # ax.set_ylabel('Altitude (km)')
    # ax.set_xlabel('Mean Cloud  Fraction ')
    # ax.set_title ('Cloud Fraction vs Altitude')
    # plt.grid(True)
    # plt.show()



    stat = 'mean'

    cl_t_g, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat.flatten(), cl_alt_lat.flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))
    cl_t_so, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat[:,constants.so_idx_1:constants.so_idx_2].flatten(), cl_alt_lat[:,constants.so_idx_1:constants.so_idx_2].flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))

    #-------------create liquid and ice cloud fraction datasets---------------#

    if 'CAM6' in directory or 'CM4' in directory or 'CM6A' in directory:
        with Dataset( constants.variable_to_filename( 'clcalipsoliq' ), 'r') as f:
            clw = np.nanmean( constants.extract_data_over_time( 'clcalipsoliq', f, start, end ), axis = 0 ) / 100 # average over time
        clw[ clw > 1 ] = np.nan

        # clw = clw[5:]
        clw_alt_lat = np.nanmean( clw, axis = -1 ) 

        with Dataset( constants.variable_to_filename( 'clcalipsoice' ), 'r') as f:
            cli = np.nanmean( constants.extract_data_over_time( 'clcalipsoice', f, start, end ), axis = 0 ) / 100 # average over time
        cli[ cli > 1 ] = np.nan
        # cli = cli[5:]
 
        cli_alt_lat = np.nanmean( cli, axis = -1 ) 

        if raw_lat.shape[0] == constants.lat.shape[0]:
            pass
        else:
            if 'CM6A' in directory:
                interpolated = interpolate.interp2d( raw_lat[5:], raw_alt, clw_alt_lat[:,5:], kind = 'cubic')
                clw_alt_lat = interpolated( constants.lat, raw_alt )

                interpolated = interpolate.interp2d( raw_lat[5:], raw_alt, cli_alt_lat[:,5:], kind = 'cubic')
                cli_alt_lat = interpolated( constants.lat, raw_alt )
            else:
                interpolated = interpolate.interp2d( raw_lat, raw_alt, clw_alt_lat, kind = 'cubic')
                clw_alt_lat = interpolated( constants.lat, raw_alt )

                interpolated = interpolate.interp2d( raw_lat, raw_alt, cli_alt_lat, kind = 'cubic')
                cli_alt_lat = interpolated( constants.lat, raw_alt )

        imp = SimpleImputer(missing_values=np.nan, strategy='mean')
        imp.fit(np.transpose(clw_alt_lat))  
        a = imp.transform(np.transpose(clw_alt_lat))
        clw_alt_lat = np.transpose(a)    

        imp = SimpleImputer(missing_values=np.nan, strategy='mean')
        imp.fit(np.transpose(cli_alt_lat))  
        a = imp.transform(np.transpose(cli_alt_lat))
        cli_alt_lat = np.transpose(a)    

        clw[ clw < 1 ] = np.nan
        cli[ cli < 1 ] = np.nan

        ######## Binned Temperature Data ########
        # values to bin: clw_frac_alt_lat and ta_alt_lat
        # binned into constants.ta_g array
        # values in each bin to be summed
        # call the summed values clw_frac_t_g

        clw_frac_t_g, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat.flatten(), clw_alt_lat.flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))
        clw_frac_t_so, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat[:,constants.so_idx_1:constants.so_idx_2].flatten(), clw_alt_lat[:,constants.so_idx_1:constants.so_idx_2].flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))

        cli_frac_t_g, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat.flatten(), cli_alt_lat.flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))
        cli_frac_t_so, bin_edges, binnumber = stats.binned_statistic(full_ta_alt_lat[:,constants.so_idx_1:constants.so_idx_2].flatten(), cli_alt_lat[:,constants.so_idx_1:constants.so_idx_2].flatten(), stat, bins=constants.ta.size, range=(constants.min_ta, constants.max_ta))


        clw_frac_g = constants.globalalt_latMean(clw_alt_lat, constants.lat)
        clw_frac_so = constants.globalalt_latMean(clw_alt_lat[:,start_idx:end_idx], constants.lat[start_idx:end_idx])
        clw_frac_g = clw_frac_g[:constants.liq_alt_confine]
        clw_frac_so = clw_frac_so[:constants.liq_alt_confine]

        cli_frac_g = constants.globalalt_latMean(cli_alt_lat, constants.lat)
        cli_frac_so = constants.globalalt_latMean(cli_alt_lat[:,start_idx:end_idx], constants.lat[start_idx:end_idx])
        

        #----Test Plot----#

        # fig, ax = plt.subplots()
        # ax.plot( clw_frac_so, constants.liq_alt )
        # ax.plot( clw_frac_g, constants.liq_alt )
        # ax.set_ylabel('Altitude (km)')
        # ax.set_xlabel('Mean Cloud  Fraction ')
        # ax.set_title ('Cloud Fraction vs Altitude')
        # plt.grid(True)
        # plt.show()



        full_clw_frac_alt_lat = clw_alt_lat
        cli_frac_alt_lat = cli_alt_lat
        clw_frac_alt_lat = full_clw_frac_alt_lat[:constants.liq_alt_confine]


        #----Test Plots----#
        # fig, ax = plt.subplots()
        # cont = ax.contourf( constants.lat, constants.liq_alt, clw_frac_alt_lat )
        # # temp = ax.contour( constants.lat[constants.lat_confine_1:constants.lat_confine_2], constants.alt, (full_ta_alt_lat[constants.lat_confine_1:constants.lat_confine_2] - 273.15), colors='white')
        # # temp.collections[5].set_linewidth(3)
        # # temp.collections[5].set_color('white')
        # # ax.clabel(temp, inline=1, fontsize=10)
        # ax.set_xlabel('Latitude')
        # ax.set_ylabel('Altitude (km)')
        # cbar = fig.colorbar(cont, orientation='horizontal')
        # cbar.set_label('Mean Cloud Liquid Water Mass Fraction in Air')
        # plt.show()


    # fig, ax = plt.subplots()
    # ax.plot( constants.ta, cl_t_g )
    # ax.plot( constants.ta, cl_t_so )
    # ax.axvline(x=273, label = '273K', color = 'black', linestyle='--')
    # plt.grid(True)
    # plt.show()

    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    ax.coastlines()
    p = ax.contourf(constants.lon, constants.lat, clt_lat_lon, transform=ccrs.PlateCarree(), cmap='coolwarm')
    cbar = plt.colorbar(p, orientation='horizontal')
    cbar.set_label('Cloud Fraction')
    ax.set_title('Total Cloud Fraction')
    plt.show()  


    ######################


    os.chdir(save_location)

    save_filename = 'COSP_' + start.strftime( '%b_%Y') + '_' + end.strftime( '%b_%Y') + '_' + directory + '.h5'

    with h5py.File( save_filename, 'w' ) as p:
        
        p.create_dataset('clt', data=clt) # total cloud fraction corresponding to lat
        p.create_dataset('clt_l', data=clt_l) # low cloud fraction corresponding to lat
        p.create_dataset('clt_m', data=clt_m) # medium cloud fraction corresponding to lat
        p.create_dataset('clt_h', data=clt_h) # high cloud fraction corresponding to lat
        p.create_dataset('clt_lat_lon', data=clt_lat_lon) # total cloud fraction corresponding to lat, lon
        p.create_dataset('clt_l_lat_lon', data=clt_l_lat_lon) # low cloud fraction corresponding to lat, lon
        p.create_dataset('clt_m_lat_lon', data=clt_m_lat_lon) # medium cloud fraction corresponding to lat, lon
        p.create_dataset('clt_h_lat_lon', data=clt_h_lat_lon) # high cloud fraction corresponding to lat, lon
      
        p.create_dataset('cl_g', data=cl_g) # global layer total cloud fraction corresponding to alt
        p.create_dataset('cl_so', data=cl_so) # southern ocean layer total cloud fraction corresponding to alt
        p.create_dataset('cl_t_g', data=cl_t_g) # global layer cloud fraction corresponding to ta
        p.create_dataset('cl_t_so', data=cl_t_so) # global layer cloud fraction corresponding to ta
        p.create_dataset('ta_alt_lat', data=ta_alt_lat) # temperature (K) corresponding to liq_alt confined to 7km and lat
        p.create_dataset('full_ta_alt_lat', data=full_ta_alt_lat) # temperature (K) corresponding to alt and lat
        p.create_dataset('cl_alt_lat', data=cl_alt_lat) # cloud fraction (%) corresponding to alt and lat

        if 'CAM6' in directory or 'CM4' in directory or 'CM6A' in directory:

            p.create_dataset('clw_frac_g', data=clw_frac_g) # global layer cloud liquid water fraction(kg/kg) corresponding to liq_alt
            p.create_dataset('cli_frac_g', data=cli_frac_g) # global layer cloud ice water fraction(kg/kg) corresponding to alt
        
            p.create_dataset('clw_frac_so', data=clw_frac_so) # southern ocean layer cloud liquid water fraction(%) corresponding to liq_alt
            p.create_dataset('cli_frac_so', data=cli_frac_so) # southern ocean layer cloud ice water fraction(%) corresponding to alt
    
            p.create_dataset('clw_frac_t_g', data=clw_frac_t_g) # global layer cloud liquid water fraction(%) corresponding to ta
            p.create_dataset('clw_frac_t_so', data=clw_frac_t_so) # global layer cloud liquid water fraction(%) corresponding to ta
    
            p.create_dataset('cli_frac_t_g', data=cli_frac_t_g) # global layer cloud ice water fraction(%) corresponding to ta
            p.create_dataset('cli_frac_t_so', data=cli_frac_t_so) # global layer cloud ice water fraction(%) corresponding to ta

            p.create_dataset('full_clw_frac_alt_lat', data=full_clw_frac_alt_lat) # cloud liquid water fraction (%) corresponding to alt and lat
            p.create_dataset('clw_frac_alt_lat', data=clw_frac_alt_lat) # cloud liquid water fraction (%) corresponding to liq_alt and lat
            p.create_dataset('cli_frac_alt_lat', data=cli_frac_alt_lat) # cloud ice liquid water fraction (%) corresponding to alt and lat

        p.close()
