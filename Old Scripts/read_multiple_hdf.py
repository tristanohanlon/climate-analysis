# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: toha006
"""

import datetime as dt  # Python standard library datetime  module
import numpy as np
import os
import tables
from pyhdf import SD
import pprint

os.chdir('Data/CERES/2016/') # change the working directory (where your files are stored)

for filename in os.listdir():
    print (filename)
    #
    f = SD.SD(filename)

#f = SD.SD('Data/CERES/CER_SSF1deg-Month_Terra-MODIS_Edition4A_401405.201601.hdf')
#view datasets
#for filename in file_names:
#    datasets_dic = f.datasets()
#    for idx,sds in enumerate(datasets_dic.keys()):
 #       print (idx,sds)
    
#sds_obj = f.select('cld_phase_37um_day_reg') # select sds

#data = sds_obj.get() # get sds data
#print (data)
#pprint.pprint( sds_obj.attributes() )