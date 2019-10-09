# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 07:25:36 2019

@author: Tristan


    CMIP5-CESM1-CAM5 : 1979.01-2005.12 --- set cl_is_fractional=False
    CMIP5-GFDL-HIRAM-C360 : 1999.01-2008.12
    CMIP5-GISS-E2R : 1951.01-2010.12 
    CMIP5-IPSL-CM5A-LR : 1979.01-2009.12
    CMIP5-MIROC5 : 1999.01-2008.12
    CMIP5-MRI-CGCM3 : 1999.01-2010.02
    
    CMIP6-CESM2-CAM6 : 1950.01-2014.12   
    CMIP6-GFDL-AM4 : 1980.01-2014.12
    CMIP6-GISS-E21G : 2001.01-2014.12
    CMIP6-IPSL-CM6A-LR : 1979.01-2014.12
    CMIP6-MIROC6 : 1999.01-2014.12
    CMIP6-MRI-ESM2 : 1999.01-2014.12

reduce_dataset is called by:
    (data folder, file name after quantity -referenced in dictionary below, destination directory - specified by location, start date, end date)

"""
import os
import datetime
import reduce_dataset
import constants

#--- Set Location and model ---#

model = 'CMIP6-CESM2-CAM6' # see above
location = constants.home # home, uni, hdd or laptop
os.chdir( location + 'Data' )

reduce_dataset.reduce_dataset( model, constants.model_dict[ model ], location + 'climate-analysis/reduced_data', datetime.datetime( 2006, 1, 1 ), datetime.datetime( 2010, 12, 1 ) )
