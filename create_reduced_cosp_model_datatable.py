# -*- coding: utf-8 -*-
"""
Reduces all data into global average data to be stored and analysed on an Excel sheet.

Calls reduce_cosp_datatable

    CMIP5-AMIP-CESM1-CAM5 : 1979.01-2005.12 --- cl must be dived by 100
    CMIP5-AMIP-GFDL-CM3 : 1999.01-2008.12
    CMIP5-AMIP-GISS-E2R : 1951.01-2010.12 
    CMIP5-AMIP-IPSL-CM5A-LR : 1979.01-2009.12
    CMIP5-AMIP-MIROC5 : 1999.01-2008.12
    CMIP5-AMIP-MRI-CGCM3 : 1999.01-2010.02
    
    CMIP6-AMIP-CESM2-CAM6 : 1950.01-2014.12   
    CMIP6-AMIP-GFDL-CM4 : 1980.01-2014.12
    CMIP6-AMIP-IPSL-CM6A-LR : 1979.01-2014.12
    CMIP6-AMIP-MIROC6 : 1999.01-2014.12
    CMIP6-AMIP-MRI-ESM2 : 1999.01-2014.12

"""
import os
import datetime
import reduce_cosp_datatable
import constants
from pprint import pprint

#--- Set Location and model ---#

start_cmip5 = datetime.datetime( 2002, 1, 1 )
end_cmip5 = datetime.datetime( 2005, 12, 1 )

start_cmip6 = datetime.datetime( 2007, 1, 1 )
end_cmip6 = datetime.datetime( 2010, 12, 1 )

#--- Set Location and model ---#

model = 'all_models' # see comments above for available models
location = constants.home # home, uni, hdd or laptop
os.chdir( location + 'Data' )

rownum = 3
if model == 'all_models':
    for model in constants.all_cosp_models:
        os.chdir( location + 'Data' )
        print('Currently reducing: ' + model)
        if 'CMIP5' in model: 
            reduce_cosp_datatable.reduce_cosp_datatable( model, location + 'climate-analysis/reduced_data', start_cmip5, end_cmip5, rownum, location )
        else:
            reduce_cosp_datatable.reduce_cosp_datatable( model, location + 'climate-analysis/reduced_data', start_cmip6, end_cmip6, rownum, location )
        rownum += 1

