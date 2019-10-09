# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon


CCCM data_types:
    'Cloud fraction profile',
    'Liquid water content profile used'
    'Ice water content profile used'
    'Temperature profile'
    
    'Irradiance layer center height profile'
    'Layer center height profile (clouds and aerosol)' 
    'Irradiance level height profile'
    'Colatitude of subsatellite point at surface at observation'
    
MISR data_types:
    'CorrCloudTopHeightFraction_Avg'
    
CERES data_types:
    'latitude'
    'cld_eff_hgt_glob'
    'cld_amount_zon'
    'cld_amount_glob'
    'cld_lwp_zon'
    'cld_iwp_zon'


"""

import numpy as np
from pyhdf import SD
import pprint
import constants

#specify location and data source - stored in constants
location = constants.hdd + '/Data/'
data = 'MISR'
data_type = 'CorrCloudTopHeightFraction_Avg'

f = SD.SD( location + constants.satellite_dict[ data ])

#view datasets
datasets_dic = f.datasets()
for idx,sds in enumerate(datasets_dic.keys()):
    print (idx,sds)
    
sds_obj = f.select( data_type ) # select sds
data = sds_obj.get() # get sds data
print (data.shape) # print data dimensions

pprint.pprint( sds_obj.attributes()['_FillValue'] ) # print data attributes
