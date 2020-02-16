# -*- coding: utf-8 -*-
"""
Reduces all data into global average data to be stored and analysed on an Excel sheet.

Calls reduce_datatable

"""
import os
import datetime
import reduce_datatable
import constants
from pprint import pprint

#--- Set Location and model ---#

start_cmip5 = datetime.datetime( 2002, 1, 1 )
end_cmip5 = datetime.datetime( 2005, 12, 1 )

start_cmip5_rcp = datetime.datetime( 2096, 1, 1 )
end_cmip5_rcp = datetime.datetime( 2099, 12, 1 )

start_cmip6 = datetime.datetime( 2007, 1, 1 )
end_cmip6 = datetime.datetime( 2010, 12, 1 )

start_cmip6_ssp = datetime.datetime( 2096, 1, 1 )
end_cmip6_ssp = datetime.datetime( 2099, 12, 1 )

#--- Set Location and model ---#

model = 'all_models' # see comments above for available models
location = constants.home # home, uni, hdd or laptop
os.chdir( location + 'Data' )

rownum = 3
if model == 'all_models':
    for model in constants.all_models:
        os.chdir( location + 'Data' )
        print('Currently reducing: ' + model)
        if 'CMIP5-AMIP' in model: 
            reduce_dataset.reduce_dataset( model, constants.model_dict_all[ model ], location + 'climate-analysis/reduced_data', start_cmip5, end_cmip5, rownum )
        elif 'CMIP5-RCP45' in model:
            reduce_dataset.reduce_dataset( model, constants.model_dict_all[ model ], location + 'climate-analysis/reduced_data', start_cmip5_rcp, end_cmip5_rcp, rownum )
        elif 'CMIP6-SSP245' in model:
            reduce_dataset.reduce_dataset( model, constants.model_dict_all[ model ], location + 'climate-analysis/reduced_data', start_cmip6_ssp, end_cmip6_ssp, rownum )
        else:
            reduce_dataset.reduce_dataset( model, constants.model_dict_all[ model ], location + 'climate-analysis/reduced_data', start_cmip6, end_cmip6, rownum )
        rownum += 1
else:
        if 'CMIP5-AMIP' in model: 
            reduce_dataset.reduce_dataset( model, constants.model_dict_all[ model ], location + 'climate-analysis/reduced_data', start_cmip5, end_cmip5, rownum )
        elif 'CMIP5-RCP45' in model:
            reduce_dataset.reduce_dataset( model, constants.model_dict_all[ model ], location + 'climate-analysis/reduced_data', start_cmip5_rcp, end_cmip5_rcp, rownum )
        elif 'CMIP6-SSP245' in model:
            reduce_dataset.reduce_dataset( model, constants.model_dict_all[ model ], location + 'climate-analysis/reduced_data', start_cmip6_ssp, end_cmip6_ssp, rownum )
        else:
            reduce_dataset.reduce_dataset( model, constants.model_dict_all[ model ], location + 'climate-analysis/reduced_data', start_cmip6, end_cmip6, rownum )

