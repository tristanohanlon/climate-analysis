# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

    CMIP5-CESM1-CAM5
    CMIP5-GFDL-HIRAM-C360
    CMIP5-GISS-E2R
    CMIP5-IPSL-CM5A-LR
    CMIP5-MIROC5
    CMIP5-MRI-CGCM3
    
    CMIP6-CESM2-CAM6
    CMIP6-GFDL-AM4
    CMIP6-GISS-E21G
    CMIP6-IPSL-CM6A-LR
    CMIP6-MIROC6
    CMIP6-MRI-ESM2
    
    ECMWF
    
    CALIPSO-GOCCP:
        3D_CloudFraction_OPAQ330m
        3D_CloudFraction_Phase330m
        3D_CloudFraction_Temp330m
        3D_CloudFraction330m
        Map_OPAQ330m
        MapLowMidHigh_Phase330m
        MapLowMidHigh330m
"""

import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import os
import constants
import numpy as np
import pprint

#specify location, data source - stored in constants and data type (cl, clw, cli, ps, ta ...)
data = 'CMIP5-RCP45-GFDL-CM3'
location = constants.home+ 'Data/'
data_type = 'clt'

if data == 'ECMWF':
    with Dataset(location + data + '/' + constants.model_dict[ data ], 'r') as f: 
        print(f.variables.keys())
        data = f.variables['level'][:]
        print(data.shape)

else:   
    with Dataset(location + data + '/' + data_type + constants.model_dict[ data ], 'r') as f: 
        print(f.variables.keys())
        data = f.variables[ 'clt' ][:]
        print(data.shape)
        time = f.variables['time'][:]
        print(time)

        # print(date2index(datetime.datetime(2096,1,1), time, select='before'))


