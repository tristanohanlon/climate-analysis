# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon
"""

import datetime as dt  # Python standard library datetime  module
import numpy as np
from pyhdf import SD
import pprint

# Specify hdf file
f = SD.SD('CCCM/2010/CER-NEWS_CCCM_Aqua-FM3-MODIS-CAL-CS_RelB1_905905.20100201.hdf')

#view datasets
datasets_dic = f.datasets()
for idx,sds in enumerate(datasets_dic.keys()):
    print (idx,sds)
    
sds_obj = f.select('Mean CloudSat radar only ice water content') # select sds

data = sds_obj.get() # get sds data
print (data.shape) # print data dimensions

pprint.pprint( sds_obj.attributes() ) # print data attributes