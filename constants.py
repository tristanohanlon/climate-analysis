# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 10:06:47 2019

@author: Tristan O'Hanlon - University of Auckland
"""
import numpy as np
from scipy import ndimage as nd

model_dict = {
    "CMIP5-CESM1-CAM5" : "_Amon_CESM1-CAM5_amip_r1i1p1_197901-200512.nc",
    "CMIP5-GFDL-HIRAM-C360" : "_Amon_GFDL-HIRAM-C360_amip_r1i1p1_199901-200812.nc",
    "CMIP5-GISS-E2R" : "_Amon_GISS-E2-R_amip_r1i1p1_195101-201012.nc",
    "CMIP5-IPSL-CM5A-LR" : "_Amon_IPSL-CM5A-LR_amip_r1i1p1_197901-200912.nc",
    "CMIP5-MIROC5" : "_Amon_MIROC5_amip_r1i1p1_199901-200812.nc",
    "CMIP5-MRI-CGCM3" : "_Amon_MRI-CGCM3_amip_r1i1p1_199901-201002.nc",
    
    "CMIP6-CESM2-CAM6" : "_Amon_CESM2_amip_r1i1p1f1_gn_195001-201412.nc",
    "CMIP6-GFDL-AM4" : "_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc",
    "CMIP6-GISS-E21G" : "_Amon_GISS-E2-1-G_amip_r1i1p1f1_gn_200101-201412.nc",
    "CMIP6-IPSL-CM6A-LR" : "_Amon_IPSL-CM6A-LR_amip_r1i1p1f1_gr_197901-201412.nc",
    "CMIP6-MIROC6" : "_Amon_MIROC6_amip_r1i1p1f1_gn_199901-201412.nc",
    "CMIP6-MRI-ESM2" : "_Amon_MRI-ESM2-0_amip_r1i1p1f1_gn_199901-201412.nc",
    
    "ECMWF" : "pressure_levels/2008_ECMWF_amon_plevels_T_cc_clw_ciw.nc",
    "CALIPSO-GOCCP" : "_200606-201803_avg_CFMIP2_sat_3.1.2.nc"
   
        }

satellite_dict = {
    "CCCM" : "CCCM/CER-NEWS_CCCM_Aqua-FM3-MODIS-CAL-CS_RelB1_905905.20060731.hdf",
    "MISR" : "MISR/MIL3YCFA/MISR_AM1_CFbA_2000.hdf",
    "CERES" : "CERES/Run/CER_SSF1deg-Month_Terra-MODIS_Edition4A_400405.200707"
        }


home = 'E:/University/University/MSc/Models/'
uni = 'C:/Users/toha006/University/University/MSc/Models/'
hdd = 'D:/MSc/Models/'
laptop = 'C:/Users/tristan/University/University/MSc/Models/'


min_lat = -75.0
max_lat = 75.0
lat_division = 0.5
lat = np.arange(min_lat, max_lat, lat_division)

# if reducing CALIPSO data, these need to be changed to min = -180, max = 180
min_lon = 0
max_lon = 360
lon_division = 1
lon = np.arange(min_lon, max_lon, lon_division)

min_alt = 0.5
max_alt = 20
alt_division = 0.5
alt = np.arange(min_alt, max_alt, alt_division)

min_ta_g = 250
max_ta_g = 280
ta_g_division = 0.5
ta_g = np.arange(min_ta_g, max_ta_g, ta_g_division)

min_ta_so = 240
max_ta_so = 270
ta_so_division = 0.5
ta_so = np.arange(min_ta_so, max_ta_so, ta_so_division)

min_liq_alt = 0.5
max_liq_alt = 7
liq_alt_division = 0.25
liq_alt = np.arange(min_liq_alt, max_liq_alt, liq_alt_division)


def fit_2d_liq_data ( data, raw_lat, raw_alt ):
    new_liq_data = np.zeros(( lat.size, liq_alt.size ))
    liq_data_counts = np.zeros( (lat.size, liq_alt.size ))
    
    for l_index, l in enumerate( raw_lat ):
        if l <= min_lat or l >= max_lat:
            continue
        lat_bin = int( ( l - min_lat ) / lat_division)
        for a_index, a in enumerate( raw_alt ):
            if a < min_liq_alt or a > max_liq_alt:
                continue
            liq_alt_bin = int( ( a - min_liq_alt ) / liq_alt_division )
            val = data[l_index, a_index]
            if np.isnan(val):
                continue
            new_liq_data[ lat_bin, liq_alt_bin ] += val
            liq_data_counts[ lat_bin, liq_alt_bin ] += 1
    data = new_liq_data / liq_data_counts
    return data
          
  
def fit_2d_data( data, raw_lat, raw_alt ):
    new_data = np.zeros(( lat.size, alt.size ))
    data_counts = np.zeros( (lat.size, alt.size ))
    
    for l_index, l in enumerate( raw_lat ):
        if l <= min_lat or l >= max_lat:
            continue
        lat_bin = int( ( l - min_lat ) / lat_division)
        for a_index, a in enumerate( raw_alt ):
            if a < min_alt or a > max_alt:
                continue
            alt_bin = int( ( a - min_alt ) / alt_division )
            val = data[l_index, a_index]
            if np.isnan(val):
                continue
            new_data[ lat_bin, alt_bin ] += val
            data_counts[ lat_bin, alt_bin ] += 1
    data = new_data / data_counts
    return data
            
           

def fit_data( data, raw_alt ):
    new_data = np.zeros(( alt.size ))
    data_counts = np.zeros( ( alt.size ))
    for a_index, a in enumerate( raw_alt ):
        if a < min_alt or a > max_alt:
            continue
        alt_bin = int( ( a - min_alt ) / alt_division )
        val = data[ a_index ]
        if np.isnan( val ):
            continue
        new_data[ alt_bin ] += val
        data_counts[ alt_bin ] += 1
            
    data = new_data / data_counts
    return data


def fit_liq_data( data, raw_alt ):
    new_data = np.zeros(( liq_alt.size ))
    data_counts = np.zeros( ( liq_alt.size ))
    for a_index, a in enumerate( raw_alt ):
        if a < min_liq_alt or a > max_liq_alt:
            continue
        liq_alt_bin = int( ( a - min_liq_alt ) / liq_alt_division )
        val = data[ a_index ]
        if np.isnan( val ):
            continue
        new_data[ liq_alt_bin ] += val
        data_counts[ liq_alt_bin ] += 1
            
    data = new_data / data_counts
    return data


def fit_lat_data( data, raw_lat ):
    new_data = np.zeros(( lat.size ))
    data_counts = np.zeros( (lat.size ))
    for l_index, l in enumerate( raw_lat ):
        if l <= min_lat or l >= max_lat:
            continue
        lat_bin = int( ( l - min_lat ) / lat_division)
        val = data[l_index]
        if np.isnan(val):
            continue
        new_data[ lat_bin ] += val
        data_counts[ lat_bin ] += 1
            
    data = new_data / data_counts
    return data

def fit_ta_g_data( data, raw_ta ):
    new_data = np.zeros(( ta_g.size ))
    data_counts = np.zeros( (ta_g.size ))
    for l_index, l in enumerate( raw_ta ):
        if l <= min_ta_g or l >= max_ta_g:
            continue
        ta_g_bin = int( ( l - min_ta_g ) / ta_g_division)
        val = data[l_index]
        if np.isnan(val):
            continue
        new_data[ ta_g_bin ] += val
        data_counts[ ta_g_bin ] += 1
            
    data = new_data / data_counts
    return data

def fit_ta_so_data( data, raw_ta ):
    new_data = np.zeros(( ta_so.size ))
    data_counts = np.zeros( (ta_so.size ))
    for l_index, l in enumerate( raw_ta ):
        if l <= min_ta_so or l >= max_ta_so:
            continue
        ta_so_bin = int( ( l - min_ta_so ) / ta_so_division)
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

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

