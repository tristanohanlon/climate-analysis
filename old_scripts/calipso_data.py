# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 18:40:37 2019

@author: Tristan O'Hanlon

This will create a dataset of time averaged global total cloud cover with latitude.
Time period is from 01.1980 - 12.2014

"""
import time
import sys
import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import h5py
import math
from scipy import integrate
import constants
from scipy import interpolate
from scipy import ndimage as nd



###########################---get latitude, longitude and cloud fraction data---###########################

os.chdir('E:/University/University/MSc/Models/Data/CALIPSO-GOCCP') #Home PC


def fit_full_data( data ):
    new_data = np.zeros(( constants.lat.size, constants.alt.size ))
    data_counts = np.zeros( (constants.lat.size, constants.alt.size ))
    for l_index, l in enumerate( raw_lat ):
        if l <= constants.min_lat or l >= constants.max_lat:
            continue
        lat_bin = int( ( l - constants.min_lat ) / constants.lat_division)
        for a_index, a in enumerate( raw_alt ):
            if a < constants.min_alt or a > constants.max_alt:
                continue
            alt_bin = int( ( a - constants.min_alt ) / constants.alt_division )
            val = data[l_index, a_index]
            if np.isnan(val):
                continue
            new_data[ lat_bin, alt_bin ] += val
            data_counts[ lat_bin, alt_bin ] += 1
            
    data = new_data / data_counts
    return data

def fit_liq_data( data ):
    new_data = np.zeros(( constants.lat.size, constants.liq_alt.size ))
    data_counts = np.zeros( (constants.lat.size, constants.liq_alt.size ))
    for l_index, l in enumerate( raw_lat ):
        if l <= constants.min_lat or l >= constants.max_lat:
            continue
        lat_bin = int( ( l - constants.min_lat ) / constants.lat_division)
        for a_index, a in enumerate( raw_alt ):
            if a < constants.min_liq_alt or a > constants.max_liq_alt:
                continue
            alt_bin = int( ( a - constants.min_liq_alt ) / constants.liq_alt_division )
            val = data[l_index, a_index]
            if np.isnan(val):
                continue
            new_data[ lat_bin, alt_bin ] += val
            data_counts[ lat_bin, alt_bin ] += 1
            
    data = new_data / data_counts
    return data

def fit_lat_data( data ):
    new_data = np.zeros(( constants.lat.size ))
    data_counts = np.zeros( (constants.lat.size ))
    for l_index, l in enumerate( raw_lat ):
        if l <= constants.min_lat or l >= constants.max_lat:
            continue
        lat_bin = int( ( l - constants.min_lat ) / constants.lat_division)
        val = data[l_index]
        if np.isnan(val):
            continue
        new_data[ lat_bin ] += val
        data_counts[ lat_bin ] += 1
            
    data = new_data / data_counts
    return data

def fit_ta_g_data( data ):
    new_data = np.zeros(( constants.ta_g.size ))
    data_counts = np.zeros( (constants.ta_g.size ))
    for l_index, l in enumerate( raw_ta ):
        if l <= constants.min_ta_g or l >= constants.max_ta_g:
            continue
        ta_g_bin = int( ( l - constants.min_ta_g ) / constants.ta_g_division)
        val = data[l_index]
        if np.isnan(val):
            continue
        new_data[ ta_g_bin ] += val
        data_counts[ ta_g_bin ] += 1
            
    data = new_data / data_counts
    return data

def fit_ta_so_data( data ):
    new_data = np.zeros(( constants.ta_so.size ))
    data_counts = np.zeros( (constants.ta_so.size ))
    for l_index, l in enumerate( raw_ta ):
        if l <= constants.min_ta_so or l >= constants.max_ta_so:
            continue
        ta_so_bin = int( ( l - constants.min_ta_so ) / constants.ta_so_division)
        val = data[l_index]
        if np.isnan(val):
            continue
        new_data[ ta_so_bin ] += val
        data_counts[ ta_so_bin ] += 1
            
    data = new_data / data_counts
    return data

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


f = Dataset('Map_OPAQ330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')
raw_lat = np.array(f.variables['latitude'][:])

#opaque (thick) clouds
cc = np.array(f.variables['cltcalipso_opaque'][:60]) # 7.2006 to 12.2020 - 54 months
cc[cc < 0] = None #set fill values to nan

cc = np.nanmean(cc, axis = 0) #average over time

#thin clouds
tf = np.array(f.variables['cltcalipso_thin'][:60]) # 7.2006 to 12.2020 - 54 months
tf[tf < 0] = None #set fill values to nan
tf = np.nanmean(tf, axis = 0) #average over time

#total clouds
cc = cc + tf #(lat, lon)
clt = np.nanmean(cc, axis = -1) #average over longitude (lat)


###########################---get alt - cloud fraction---###########################

f = Dataset('3D_CloudFraction330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')
raw_alt = np.array(f.variables['alt_mid'][:]) # km

cl = np.array(f.variables['clcalipso'][:60]) # 7.2006 to 12.2020 - 54 months

cl[cl < 0] = None #set fill values to nan

cl = np.nanmean(cl, axis = 0) #average over time
cl_alt_lat = np.nanmean(cl, axis = -1) #average over longitude
cl_g = np.nanmean(cl_alt_lat, axis = -1) #average over latitude


###########################---get alt - phase fractions---###########################

f = Dataset('3D_CloudFraction_Phase330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')

#ice water fraction
ice = np.array(f.variables['clcalipso_RPIC'][0:60]) # 2.2011 to 4.2011 - 3 months
ice[ice < 0] = None #set fill values to nan
ice = np.nanmean(ice, axis = 0) #average over time

#ice water fraction with alt and lat

cli_alt_lat = np.nanmean(ice, axis = -1) #average over lon
clw_alt_lat = (1 - cli_alt_lat) * cl_alt_lat
cli_alt_lat = cli_alt_lat * cl_alt_lat 
#ice water fraction
cli_g = np.nanmean(cli_alt_lat, axis = -1) #average over lat
clw_g = np.nanmean(clw_alt_lat, axis = -1) #average over lat


###############---create southern ocean data---###############


cl_so = np.hstack((np.vstack(lat), np.transpose(cl_alt_lat) )) #creates a (180,34) array
cl_so = cl_so[cl_so[:,0]>=-70]
cl_so = cl_so[cl_so[:,0]<=-50]
cl_so = cl_so[:,1:] #Split the combined array into just the tccl data, eliminating the first coloumn of latitude
cl_so = np.nanmean(cl_so, axis = 0)


clw_so = np.hstack((np.vstack(lat), np.transpose(clw_alt_lat))) #creates a (180,34) array
clw_so = clw_so[clw_so[:,0]>=-70]
clw_so = clw_so[clw_so[:,0]<=-50]
clw_so = clw_so[:,1:] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
clw_so = np.nanmean(clw_so, axis = 0)


cli_so = np.hstack((np.vstack(lat), np.transpose(cli_alt_lat))) #creates a (180,34) array
cli_so = cli_so[cli_so[:,0]>=-70]
cli_so = cli_so[cli_so[:,0]<=-50]
cli_so = cli_so[:,1:] #Split the combined array into just the tclw data, eliminating the first coloumn of latitude
cli_so = np.nanmean(cli_so, axis = 0)


###########################---get alt - temp - phase fractions---###########################

f = Dataset('3D_CloudFraction_Temp330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r')
raw_ta = np.array(f.variables['temp_mid'][:]) + 273 # km (38)

cl_t = np.array(f.variables['cltemp'][0:60]) # 7.2006 to 12.2020 - 54 months
cl_t [cl_t  < 0] = None #set fill values to nan
cl_t  = np.nanmean(cl_t, axis = 0) #average over time
cl_t_lat =  np.nanmean(cl_t, axis = -1) #average over longitude

cli_t = np.array(f.variables['cltemp_phase'][0:60]) # 7.2006 to 12.2020 - 54 months
cli_t[cli_t < 0] = None #set fill values to nan

cli_t = np.nanmean(cli_t, axis = 0) #average over time
cli_t_lat =  np.nanmean(cli_t, axis = -1) #average over longitude

clw_t_lat =  (1 - cli_t_lat) * cl_t_lat   
clw_t_g = np.nanmean(clw_t_lat, axis = -1)

clw_t_so = np.hstack((np.vstack(lat), np.transpose(clw_t_lat))) #creates a (180,34) array
clw_t_so = clw_t_so[clw_t_so[:,0]>=-70]
clw_t_so = clw_t_so[clw_t_so[:,0]<=-50]
clw_t_so = clw_t_so[:,1:] #Split the combined array into just the tcclw data, eliminating the first coloumn of latitude
clw_t_so = np.nanmean(clw_t_so, axis = 0)


clt = fit_lat_data(clt)
cl_alt_lat = fit_full_data( np.transpose(cl_alt_lat) )
full_clw_alt_lat = fit_full_data( np.transpose(clw_alt_lat) )
cli_alt_lat = fit_full_data( np.transpose(cli_alt_lat) )
clw_alt_lat = fit_liq_data( np.transpose(clw_alt_lat) )

clw_t_g = fit_ta_g_data(clw_t_g)
clw_t_so = fit_ta_so_data(clw_t_so)

###############---create combined data---###############

#os.chdir('E:/University/University/MSc/Models/climate-analysis/reduced_datasets') #Home PC
#
#with h5py.File('06.2006_06.2011_CALIPSO.h5', 'w') as p:
#    
#    p.create_dataset('lat', data=lat)  
#    p.create_dataset('alt', data=alt)
#    p.create_dataset('alt_t', data=alt_t)
#    
#    p.create_dataset('tcc', data=tcc)
#    p.create_dataset('tclw_frac', data=tclw)
#    p.create_dataset('tciw_frac', data=tciw)
# 
#    p.create_dataset('cf_t_lat', data=cf_t_lat)
#    p.create_dataset('lw_t_lat', data=lw_t_lat)
#    p.create_dataset('iw_t_lat', data=iw_t_lat)
#    
#    p.create_dataset('cf', data=cf)
#    p.create_dataset('lw_frac', data=lw_frac)
#    p.create_dataset('iw_frac', data=iw_frac)
#    
#    p.create_dataset('cf_so', data=cf_so)
#    p.create_dataset('lw_frac_so', data=lw_so)
#    p.create_dataset('iw_frac_so', data=iw_so)   
#
#    p.create_dataset('cf_t', data=cf_t)
#    p.create_dataset('lw_t_frac', data=lw_t)
#    p.create_dataset('iw_t_frac', data=iw_t)   
#
#    p.create_dataset('cf_t_so', data=cf_t_so)
#    p.create_dataset('lw_t_frac_so', data=lw_t_so)
#    p.create_dataset('iw_t_frac_so', data=iw_t_so)   
#    
#    p.create_dataset('cf_alt_lat', data=cff)
#    p.create_dataset('liq_frac_alt_lat', data=liq_frac_alt_lat)
#    p.create_dataset('ice_frac_alt_lat', data=ice_frac_alt_lat)
#    
#    p.close()



