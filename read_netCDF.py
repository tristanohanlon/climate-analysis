# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

    "CMIP5-AMIP-CESM1-CAM5" : "_cfMon_CESM1-CAM5_amip_r2i1p1_197901-200512.nc",
    "CMIP5-AMIP-GFDL-CM3" : "_cfMon_GFDL-CM3_amip_r1i1p1_199901-200812.nc",
    "CMIP5-AMIP-MIROC5" : "_cfMon_MIROC5_amip_r1i1p1_199901-200812",
    "CMIP5-AMIP-MRI-CGCM3" : "_cfMon_MRI-CGCM3_amip_r1i1p1_199901-200812.nc",

    "CMIP6-AMIP-CESM2-CAM6" : "_CFmon_CESM2_amip_r2i1p1f1_gn_195001-201412.nc",
    "CMIP6-AMIP-GFDL-CM4" : "_CFmon_GFDL-CM4_amip_r1i1p1f1_gr1_200301-201412.nc",
    "CMIP6-AMIP-MIROC6" : "_CFmon_MIROC6_amip_r1i1p1f1_gn_199901-201412.nc",
    "CMIP6-AMIP-MRI-ESM2" : "_CFmon_MRI-ESM2-0_amip_r1i1p1f1_gn_199901-201412.nc",

    
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
data = 'CMIP6-AMIP-GFDL-CM4'
location = constants.hdd + 'Data/'
data_type = 'o3'

if data == 'ECMWF':
    with Dataset(location + data + '/' + constants.model_dict[ data ], 'r') as f: 
        print(f.variables.keys())
        data = f.variables['level'][:]
        print(data.shape)

else:   
    with Dataset(location + data + data_type + constants.model_dict_cosp[ data ], 'r') as f: 
        print(f.variables.keys())
        data = f.variables[ 'alt40' ][:] / 1000
        print(data)
        # time = f.variables['time'][:]
        # print(time)

        # print(date2index(datetime.datetime(2096,1,1), time, select='before'))


