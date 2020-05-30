# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 07:25:36 2019

@author: Tristan


    CMIP5-AMIP-CESM1-CAM5 : 1979.01-2005.12 --- cl must be dived by 100
    CMIP5-AMIP-GFDL-CM3 : 1999.01-2008.12
    CMIP5-AMIP-GISS-E2R : 1951.01-2010.12 
    CMIP5-AMIP-IPSL-CM5A-LR : 1979.01-2009.12
    CMIP5-AMIP-MIROC5 : 1999.01-2008.12
    CMIP5-AMIP-MRI-CGCM3 : 1999.01-2010.02
    
    CMIP6-AMIP-CESM2-CAM6 : 1950.01-2014.12   
    CMIP6-AMIP-GFDL-CM4 : 1980.01-2014.12
    CMIP6-AMIP-GISS-E21G : 2001.01-2014.12
    CMIP6-AMIP-IPSL-CM6A-LR : 1979.01-2014.12
    CMIP6-AMIP-MIROC6 : 1999.01-2014.12
    CMIP6-AMIP-MRI-ESM2 : 1999.01-2014.12

If cosp_status is set to true, the data is taken from a secondary directory 
where the data has been run through a LIDAR simulator to be compared 
with satellite data.

reduce_dataset is called by:
    (data folder, amip file name after quantity - referenced in dictionary [model_dict_amip],
    cosp file name after quantity - referenced in dictionary [model_dict_cosp],  
        destination directory - specified by location, start date, end date, cosp status)

data is saved in a *.h5 file containing the time averages parameters - 
see the reduced dataset function for data details
"""
import os
import datetime
import reduce_dataset
import reduce_cosp_dataset
import constants
from pprint import pprint

#--- Set Location and model ---#

model = 'all_models' # see intro comments above for available models. 'all_models' runs through all
location = constants.home # home, uni, hdd or laptop
cosp_status = False # see intro comments
os.chdir( location + 'Data' ) # location of the model data folders

# set date parameters for CMIP5 and CMIP6 data

def model_date( model ):
    if 'CMIP5' in model:
        start = datetime.datetime( 2002, 1, 1 )
        end = datetime.datetime( 2005, 12, 1 )
    else:
        start = datetime.datetime( 2007, 1, 1 )
        end = datetime.datetime( 2010, 12, 1 )
    return start, end

if cosp_status == True:

    if model == 'all_models':
        for model in constants.all_cosp_models:
            (start, end) = model_date(model)
            reduce_cosp_dataset.reduce_cosp_dataset(
                model, 
                location + 'climate-analysis/reduced_data', 
                start, 
                end,
                location )

    else:
        (start, end) = model_date(model)
        reduce_cosp_dataset.reduce_cosp_dataset( 
            model, 
            location + 'climate-analysis/reduced_data', 
            start, 
            end,
            location )
else:

    if model == 'all_models':
        for model in constants.all_amip_models:
            (start, end) = model_date(model)
            reduce_dataset.reduce_dataset(
                model, 
                location + 'climate-analysis/reduced_data', 
                start, 
                end,
                location )

    else:
        (start, end) = model_date(model)
        reduce_dataset.reduce_dataset( 
            model, 
            location + 'climate-analysis/reduced_data', 
            start, 
            end,
            location )
