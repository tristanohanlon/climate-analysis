# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 10:06:47 2019

@author: Tristan O'Hanlon - University of Auckland
"""
import numpy as np
import os
from scipy import ndimage as nd
from netCDF4 import date2index

home = 'E:/University/University/MSc/Models/'
network = '//synthesis/E/University/University/MSc/Models/'
uni = 'C:/Users/toha006/University/University/MSc/Models/'
hdd = 'E:/MSc/Models/'
laptop = 'C:/Users/tristan/University/University/MSc/Models/'


all_amip_models = [ 'CMIP5-AMIP-CESM1-CAM5', 
                'CMIP5-AMIP-GFDL-CM3',
                'CMIP5-AMIP-GISS-E2R', 
                'CMIP5-AMIP-IPSL-CM5A-LR', 
                'CMIP5-AMIP-MIROC5', 
                'CMIP5-AMIP-MRI-CGCM3',

                'CMIP6-AMIP-CESM2-CAM6', 
                'CMIP6-AMIP-GFDL-CM4',
                'CMIP6-AMIP-GISS-E21G', 
                'CMIP6-AMIP-IPSL-CM6A-LR', 
                'CMIP6-AMIP-MIROC6', 
                'CMIP6-AMIP-MRI-ESM2' ]

all_cosp_models = [ 'CMIP5-AMIP-CESM1-CAM5', 
                'CMIP5-AMIP-GFDL-CM3',
                'CMIP5-AMIP-MIROC5', 
                'CMIP5-AMIP-MRI-CGCM3',

                'CMIP6-AMIP-CESM2-CAM6', 
                'CMIP6-AMIP-GFDL-CM4',
                'CMIP6-AMIP-MIROC6', 
                'CMIP6-AMIP-MRI-ESM2' ]


date_cmip5 = 'Jan_2002_Dec_2005'
date_cmip6 = 'Jan_2007_Dec_2010'
date_ceres = 'Jan_2007_Dec_2010'
date_cccm = 'Jan_2007_Dec_2010'


min_lat = -89.5
max_lat = 90
lat_division = 1
lat = np.arange(min_lat, max_lat, lat_division)
# Southern Ocean index ranges
so_idx_1 = np.abs(lat - (-70)).argmin()
so_idx_2 = np.abs(lat - (-50)).argmin()
# Limit latitude plot range to -75 and 75 degrees
lat_confine_1 = np.abs(lat - (-75)).argmin()
lat_confine_2 = np.abs(lat - (75)).argmin()

# if reducing CALIPSO data, these need to be changed to min = -180, max = 180
min_lon = 0
max_lon = 360
lon_division = 1
lon = np.arange(min_lon, max_lon, lon_division)

alt = np.array( [0.24, 0.72, 1.2, 1.68, 2.16, 2.64, 3.12, 3.6, 4.08, 4.56, 5.04, 5.52, 6.0, 6.48, 6.96, 7.44, 7.92, 8.4, 8.88, 9.36, 9.84, 10.32, 10.8, 11.28, 11.76, 12.24, 12.72, 13.2, 13.68, 14.16, 14.64, 15.12, 15.6, 16.08, 16.56, 17.04, 17.52, 18.0, 18.48, 18.96] )
liq_alt_confine = np.abs(alt - (7)).argmin()
liq_alt = alt[:liq_alt_confine]


min_ta = 240
max_ta = 300
ta_division = 2
ta = np.arange(min_ta, max_ta, ta_division)


# Function that extracts data over time.
def extract_data_over_time( type, f, start, end ):
        time_variable = f.variables['time']
        start_index = date2index( start, time_variable, select='before' )
        end_index = date2index( end, time_variable, select='before' )

        data = np.array( f.variables[type][start_index:end_index])
        return data

# Function that returns a filename based on the input variable

def variable_to_filename( variable ):
    dir_list = os.listdir() 
    for f in dir_list:
        if f.startswith(variable + '_'):
            return f



# Function that extracts all data
def extract_data( type, f ):
    return np.array( f.variables[type][:] )



# Function that selects only southern ocean data
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


# Function to calculate the global mean of a variable that has
# already been averaged over time. 
# Data is (lat)
# Input the raw latitudes from the variable file.
# Output is a single value
def globalMean(Data, latitudes): 
    areaWeights = np.cos(latitudes*np.pi/180)
    weightedMatrix = Data*areaWeights
    sumWeighted = np.sum(weightedMatrix,axis=0)
    sumWeights = np.sum(areaWeights,axis=0)
    weightedMean = sumWeighted/sumWeights
    return weightedMean


# Function to calculate the zonal low cloud mean of a variable that has
# already been averaged over time and longitude. 
# Data is (alt, lat)
# Input the raw pressures from the variable file.
# Output is (lat)
def lowregMean(Data, pressures): 
    areaWeights = pressures
    areaWeights3D = np.swapaxes(np.swapaxes(np.tile(areaWeights,
                                (np.shape(Data)[1],np.shape(Data)[2],1)), 0, 2), 1, 2)  #Replicating the area weights 
    weighted3DMatrix = Data*areaWeights3D
    sumWeighted = np.sum(weighted3DMatrix,axis=(0))
    sumWeights3D = np.sum(areaWeights3D,axis=(0))
    weightedMean = sumWeighted/sumWeights3D
    return weightedMean


def lowlatMean(Data, pressures): 
    areaWeights = pressures
    areaWeights3D = np.swapaxes(np.swapaxes(np.tile(areaWeights,
                                (np.shape(Data)[1],np.shape(Data)[2],1)), 0, 2), 1, 2)  #Replicating the area weights 
    weighted3DMatrix = Data*areaWeights3D
    sumWeighted = np.sum(weighted3DMatrix,axis=(0, 2))
    sumWeights3D = np.sum(areaWeights3D,axis=(0, 2))
    weightedMean = sumWeighted/sumWeights3D
    return weightedMean


# Function to calculate the global mean of a variable that has
# already been averaged over time. 
# Data is (lat, lon)
# Input the raw latitudes from the variable file.
# Output is a single value
def global2DMean(Data, latitudes): 
    areaWeights = np.cos(latitudes*np.pi/180)
    areaWeights2D = np.swapaxes(np.tile(areaWeights,
                                (np.shape(Data)[1], 1)), 0, 1)  #Replicating the area weights 
    weighted2DMatrix = Data*areaWeights2D
    sumWeighted = np.sum(weighted2DMatrix,axis=(0,1))
    sumWeights2D = np.sum(areaWeights2D,axis=(0,1))
    weightedMean = sumWeighted/sumWeights2D
    return weightedMean

# Function to calculate the global mean of a variable that has
# already been averaged over time. 
# Data is (alt, lat)
# Input the raw latitudes from the variable file.
# Output is a altitude array
def globalalt_latMean(Data, latitudes): 
    areaWeights = np.cos(latitudes*np.pi/180)
    areaWeights2D = np.tile(areaWeights,
                                (np.shape(Data)[0], 1))  #Replicating the area weights 
    weighted2DMatrix = Data*areaWeights2D
    sumWeighted = np.sum(weighted2DMatrix,axis=(1))
    sumWeights2D = np.sum(areaWeights2D,axis=(1))
    weightedMean = sumWeighted/sumWeights2D
    return weightedMean


# Function to calculate the global mean of a variable that has
# already been averaged over time. 
# Data is (alt, lat)
# Input the raw latitudes from the variable file.
# Output is a single value
def globalalt_latMeanVal(Data, latitudes): 
    areaWeights = np.cos(latitudes*np.pi/180)
    areaWeights2D = np.tile(areaWeights,
                                (np.shape(Data)[0], 1))  #Replicating the area weights 
    weighted2DMatrix = Data*areaWeights2D
    sumWeighted = np.sum(weighted2DMatrix,axis=(0,1))
    sumWeights2D = np.sum(areaWeights2D,axis=(0,1))
    weightedMean = sumWeighted/sumWeights2D
    return weightedMean


# Function to calculate the global profile mean of a variable that has
# already been averaged over time. 
# Data is (alt, lat, lon)
# Input the raw latitudes from the variable file.
# Output is an array which corresponds to pressure levels
def global3DMean(Data, latitudes):
    areaWeights = np.cos(latitudes*np.pi/180)
    areaWeights3D = np.swapaxes(np.tile(areaWeights,
                                (np.shape(Data)[0],np.shape(Data)[2],1)),
                                1,2)  #Replicating the area weights 
    weighted3DMatrix = Data*areaWeights3D
    sumWeighted = np.sum(weighted3DMatrix,axis=(1,2))
    sumWeights3D = np.sum(areaWeights3D,axis=(1,2))
    weightedMean = sumWeighted/sumWeights3D
    return weightedMean


# Function to calculate the global profile mean of a variable that has
# already been averaged over time. 
# Data is (alt, lat, lon)
# Input the raw latitudes from the variable file.
# Output is a single mean value
def global3DMeanVal(Data, latitudes):
    areaWeights = np.cos(latitudes*np.pi/180)
    areaWeights3D = np.swapaxes(np.tile(areaWeights,
                                (np.shape(Data)[0],np.shape(Data)[2],1)),
                                1,2)  #Replicating the area weights 
    weighted3DMatrix = Data*areaWeights3D
    sumWeighted = np.sum(weighted3DMatrix,axis=(0,1,2))
    sumWeights3D = np.sum(areaWeights3D,axis=(0,1,2))
    weightedMean = sumWeighted/sumWeights3D
    return weightedMean

# Function that bins the data into liquid altitudes and latitudes.
# Input the data that is (alt, lat) as well as the raw alts and lats
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
    

# Function that bins the data into latitudes and longitude.
# Input the data that is (lat, lon) as well as the raw lons and lats
def fit_grid_data( data, raw_lat, raw_lon ):
    new_data = np.zeros(( lat.size, lon.size ))
    data_counts = np.zeros( (lat.size, lon.size ))
    
    for l_index, l in enumerate( raw_lat ):
        if l <= min_lat or l >= max_lat:
            continue
        lat_bin = int( ( l - min_lat ) / lat_division)
        for a_index, a in enumerate( raw_lon ):
            if a < min_lon or a > max_lon:
                continue
            lon_bin = int( ( a - min_lon ) / lon_division )
            val = data[l_index, a_index]
            if np.isnan(val):
                continue
            new_data[ lat_bin, lon_bin ] += val
            data_counts[ lat_bin, lon_bin ] += 1
    data = new_data / data_counts
    return data
      
# Function that bins the data into altitudes and latitudes.
# Input the data that is (alt, lat) as well as the raw alts and lats 
def fit_2d_data( data, raw_lat, raw_alt ):
    new_data = np.zeros(( lat.size, alt.size ))
    data_counts = np.zeros( (lat.size, alt.size ))
    
    for l_index, l in enumerate( raw_lat ):
        if l <= min_lat or l >= max_lat:
            continue
        lat_bin = int( ( l - min_lat ) / lat_division)
        for a_index, a in enumerate( raw_alt ):
            if a < alt.min or a > alt.max:
                continue
            alt_bin = int( ( a - alt.min ) / alt_division )
            val = data[l_index, a_index]
            if np.isnan(val):
                continue
            new_data[ lat_bin, alt_bin ] += val
            data_counts[ lat_bin, alt_bin ] += 1
    data = new_data / data_counts
    return data
            
           
# Function that bins the data into altitudes.
# Input the data that is (alt) as well as the raw alts
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


# Function that bins the data into liquid altitudes.
# Input the data that is (alt) as well as the raw alts
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


# Function that bins the data into latitudes.
# Input the data that is (lats) as well as the raw lats
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


# Function that bins the data into temperatures.
# Input the data that is () as well as the raw ta
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

# Function that bins the data into temperatures.
# Input the data that is () as well as the raw ta
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


# Function that fills that nan vaules in the data based on nearby data
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


all_models = [ 'CMIP5-AMIP-CESM1-CAM5', 
                'CMIP5-AMIP-GFDL-HIRAM-C360', 
                'CMIP5-AMIP-GISS-E2R', 
                'CMIP5-AMIP-IPSL-CM5A-LR', 
                'CMIP5-AMIP-MIROC5', 
                'CMIP5-AMIP-MRI-CGCM3',
                'CMIP5-AMIP-GFDL-CM3',

                # 'CMIP5-AMIP_4xCO2-GFDL-CM3',
                'CMIP5-AMIP_4xCO2-IPSL-CM5A-LR',
                'CMIP5-AMIP_4xCO2-MIROC5',
                'CMIP5-AMIP_4xCO2-MRI-CGCM3',

                'CMIP5-RCP45-CESM1-BGC',
                'CMIP5-RCP45-GFDL-CM3',
                'CMIP5-RCP45-GISS-E2R',
                'CMIP5-RCP45-IPSL-CM5A-LR',
                'CMIP5-RCP45-MIROC5',
                'CMIP5-RCP45-MRI-CGCM3',

                'CMIP6-AMIP-CESM2-CAM6', 
                'CMIP6-AMIP-GFDL-AM4', 
                'CMIP6-AMIP-GISS-E21G', 
                'CMIP6-AMIP-IPSL-CM6A-LR', 
                'CMIP6-AMIP-MIROC6', 
                'CMIP6-AMIP-MRI-ESM2',
                'CMIP6-AMIP-GFDL-CM4',
                
                'CMIP6-AMIP_4xCO2-CESM2-CAM6',
                'CMIP6-AMIP_4xCO2-GFDL-CM4',
                'CMIP6-AMIP_4xCO2-IPSL-CM6A-LR',
                'CMIP6-AMIP_4xCO2-MIROC6',
                'CMIP6-AMIP_4xCO2-MRI-ESM2',

                'CMIP6-AMIP_future4K-CESM2',
                # 'CMIP6-AMIP_future4K-GFDL-CM4',
                'CMIP6-AMIP_future4K-IPSL-CM6A-LR',
                'CMIP6-AMIP_future4K-MIROC6',
                'CMIP6-AMIP_future4K-MRI-ESM2',

                'CMIP6-SSP245-CESM2',
                'CMIP6-SSP245-GFDL-CM4',
                'CMIP6-SSP245-IPSL-CM6A-LR',
                'CMIP6-SSP245-MIROC6',
                'CMIP6-SSP245-MRI-ESM2' ]


model_dict_cosp = {
    "CMIP5-AMIP-CESM1-CAM5" : "_cfMon_CESM1-CAM5_amip_r2i1p1_197901-200512.nc",
    "CMIP5-AMIP-GFDL-CM3" : "_cfMon_GFDL-CM3_amip_r1i1p1_199901-200812.nc",
    "CMIP5-AMIP-MIROC5" : "_cfMon_MIROC5_amip_r1i1p1_199901-200812.nc",
    "CMIP5-AMIP-MRI-CGCM3" : "_cfMon_MRI-CGCM3_amip_r1i1p1_199901-200812.nc",

    "CMIP6-AMIP-CESM2-CAM6" : "_CFmon_CESM2_amip_r2i1p1f1_gn_195001-201412.nc",
    "CMIP6-AMIP-GFDL-CM4" : "_CFmon_GFDL-CM4_amip_r1i1p1f1_gr1_200301-201412.nc",
    "CMIP6-AMIP-MIROC6" : "_CFmon_MIROC6_amip_r1i1p1f1_gn_199901-201412.nc",
    "CMIP6-AMIP-MRI-ESM2" : "_CFmon_MRI-ESM2-0_amip_r1i1p1f1_gn_199901-201412.nc",
}


model_dict_amip = {
    "CMIP5-AMIP-CESM1-CAM5" : "_Amon_CESM1-CAM5_amip_r1i1p1_197901-200512.nc",
    "CMIP5-AMIP-GFDL-CM3" : "_Amon_GFDL-CM3_amip_r1i1p1_199901-200812.nc",
    "CMIP5-AMIP-GISS-E2R" : "_Amon_GISS-E2-R_amip_r1i1p1_195101-201012.nc",
    "CMIP5-AMIP-IPSL-CM5A-LR" : "_Amon_IPSL-CM5A-LR_amip_r1i1p1_197901-200912.nc",
    "CMIP5-AMIP-MIROC5" : "_Amon_MIROC5_amip_r1i1p1_199901-200812.nc",
    "CMIP5-AMIP-MRI-CGCM3" : "_Amon_MRI-CGCM3_amip_r1i1p1_199901-201002.nc",

    "CMIP6-AMIP-CESM2-CAM6" : "_Amon_CESM2_amip_r1i1p1f1_gn_195001-201412.nc",
    "CMIP6-AMIP-GFDL-CM4" : "_Amon_GFDL-CM4_amip_r1i1p1f1_gr1_200301-201412.nc",
    "CMIP6-AMIP-GISS-E21G" : "_Amon_GISS-E2-1-G_amip_r1i1p1f1_gn_200101-201412.nc",
    "CMIP6-AMIP-IPSL-CM6A-LR" : "_Amon_IPSL-CM6A-LR_amip_r1i1p1f1_gr_197901-201412.nc",
    "CMIP6-AMIP-MIROC6" : "_Amon_MIROC6_amip_r1i1p1f1_gn_199901-201412.nc",
    "CMIP6-AMIP-MRI-ESM2" : "_Amon_MRI-ESM2-0_amip_r1i1p1f1_gn_199901-201412.nc",
}

model_dict_all = {
    "CMIP5-AMIP-CESM1-CAM5" : "_Amon_CESM1-CAM5_amip_r1i1p1_197901-200512.nc",
    "CMIP5-AMIP-GFDL-HIRAM-C360" : "_Amon_GFDL-HIRAM-C360_amip_r1i1p1_199901-200812.nc",
    "CMIP5-AMIP-GISS-E2R" : "_Amon_GISS-E2-R_amip_r1i1p1_195101-201012.nc",
    "CMIP5-AMIP-IPSL-CM5A-LR" : "_Amon_IPSL-CM5A-LR_amip_r1i1p1_197901-200912.nc",
    "CMIP5-AMIP-MIROC5" : "_Amon_MIROC5_amip_r1i1p1_199901-200812.nc",
    "CMIP5-AMIP-MRI-CGCM3" : "_Amon_MRI-CGCM3_amip_r1i1p1_199901-201002.nc",
    "CMIP5-AMIP-GFDL-CM3" : "_Amon_GFDL-CM3_amip_r1i1p1_199901-200812.nc",

    "CMIP5-AMIP_4xCO2-GFDL-CM3" : "_Amon_GFDL-CM3_amip4xCO2_r1i1p1_199901-200812.nc",
    "CMIP5-AMIP_4xCO2-IPSL-CM5A-LR" : "_Amon_IPSL-CM5A-LR_amip4xCO2_r1i1p1_197901-200912.nc",
    "CMIP5-AMIP_4xCO2-MIROC5" : "_Amon_MIROC5_amip4xCO2_r1i1p1_199901-200812.nc",
    "CMIP5-AMIP_4xCO2-MRI-CGCM3" : "_Amon_MRI-CGCM3_amip4xCO2_r1i1p1_199901-200812.nc",

    "CMIP5-AMIP_future4K-CESM1-CAM5" : "_Amon_CESM1-CAM5_amip4K_r1i1p1_197901-200012.nc",

    "CMIP5-RCP45-CESM1-BGC" : "_Amon_CESM1-BGC_rcp45_r1i1p1_205001-210012.nc",
    "CMIP5-RCP45-GFDL-CM3" : "_Amon_GFDL-CM3_rcp45_r1i1p1_209601-210012.nc",
    "CMIP5-RCP45-GISS-E2R" : "_Amon_GISS-E2-R_rcp45_r1i1p1_207601-210012.nc",
    "CMIP5-RCP45-IPSL-CM5A-LR" : "_Amon_IPSL-CM5A-LR_rcp45_r1i1p1_200601-210512.nc",
    "CMIP5-RCP45-MIROC5" : "_Amon_MIROC5_rcp45_r1i1p1_209001-209912.nc",
    "CMIP5-RCP45-MRI-CGCM3" : "_Amon_MRI-CGCM3_rcp45_r1i1p1_209601-210012.nc",

    "CMIP6-AMIP-CESM2-CAM6" : "_Amon_CESM2_amip_r1i1p1f1_gn_195001-201412.nc",
    "CMIP6-AMIP-GFDL-AM4" : "_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc",
    "CMIP6-AMIP-GISS-E21G" : "_Amon_GISS-E2-1-G_amip_r1i1p1f1_gn_200101-201412.nc",
    "CMIP6-AMIP-IPSL-CM6A-LR" : "_Amon_IPSL-CM6A-LR_amip_r1i1p1f1_gr_197901-201412.nc",
    "CMIP6-AMIP-MIROC6" : "_Amon_MIROC6_amip_r1i1p1f1_gn_199901-201412.nc",
    "CMIP6-AMIP-MRI-ESM2" : "_Amon_MRI-ESM2-0_amip_r1i1p1f1_gn_199901-201412.nc",
    "CMIP6-AMIP-GFDL-CM4" : "_Amon_GFDL-CM4_amip_r1i1p1f1_gr1_200301-201412.nc",

    "CMIP6-AMIP_4xCO2-CESM2-CAM6" : "_Amon_CESM2_amip-4xCO2_r1i1p1f1_gn_197901-201412.nc",
    "CMIP6-AMIP_4xCO2-GFDL-CM4" : "_Amon_GFDL-CM4_amip-4xCO2_r1i1p1f1_gr1_197901-201412.nc",
    "CMIP6-AMIP_4xCO2-IPSL-CM6A-LR" : "_Amon_IPSL-CM6A-LR_amip-4xCO2_r1i1p1f1_gr_197901-201412.nc",
    "CMIP6-AMIP_4xCO2-MIROC6" : "_Amon_MIROC6_amip-4xCO2_r1i1p1f1_gn_199901-201412.nc",
    "CMIP6-AMIP_4xCO2-MRI-ESM2" : "_Amon_MRI-ESM2-0_amip-4xCO2_r1i1p1f1_gn_199901-201412.nc",

    "CMIP6-AMIP_future4K-CESM2" : "_Amon_CESM2_amip-future4K_r1i1p1f1_gn_197901-201412.nc",
    "CMIP6-AMIP_future4K-GFDL-CM4" : "_Amon_GFDL-CM4_amip-future4K_r1i1p1f1_gr1_197901-201412.nc",
    "CMIP6-AMIP_future4K-IPSL-CM6A-LR" : "_Amon_IPSL-CM6A-LR_amip-future4K_r1i1p1f1_gr_197901-201412.nc",
    "CMIP6-AMIP_future4K-MIROC6" : "_Amon_MIROC6_amip-future4K_r1i1p1f1_gn_200901-201412.nc",
    "CMIP6-AMIP_future4K-MRI-ESM2" : "_Amon_MRI-ESM2-0_amip-future4K_r1i1p1f1_gn_200901-201412.nc",

    "CMIP6-SSP245-CESM2" : "_Amon_CESM2_ssp245_r1i1p1f1_gn_206501-210012.nc",
    "CMIP6-SSP245-GFDL-CM4" : "_Amon_GFDL-CM4_ssp245_r1i1p1f1_gr1_201501-210012.nc",
    "CMIP6-SSP245-IPSL-CM6A-LR" : "_Amon_IPSL-CM6A-LR_ssp245_r1i1p1f1_gr_201501-210012.nc",
    "CMIP6-SSP245-MIROC6" : "_Amon_MIROC6_ssp245_r1i1p1f1_gn_209501-210012.nc",
    "CMIP6-SSP245-MRI-ESM2" : "_Amon_MRI-ESM2-0_ssp245_r1i1p1f1_gn_209501-210012.nc",
    
    "ECMWF" : "pressure_levels/2008_ECMWF_amon_plevels_T_cc_clw_ciw.nc",
    "CALIPSO-GOCCP" : "_200606-201803_avg_CFMIP2_sat_3.1.2.nc"
   
        }

satellite_dict = {
    "CCCM" : "CCCM/CER-NEWS_CCCM_Aqua-FM3-MODIS-CAL-CS_RelB1_905905.20060731.hdf",
    "MISR" : "MISR/MIL3YCFA/MISR_AM1_CFbA_2000.hdf",
    "CERES" : "CERES/Run/CER_SSF1deg-Month_Terra-MODIS_Edition4A_400405.200707"
        }