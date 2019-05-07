# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon
"""

import datetime as dt  # Python standard library datetime  module
import numpy as np
from pyhdf import SD
import pprint
import os

# Specify hdf file
#f = SD.SD('C:/Users/tristan/University/University/MSc/Models/Data/CCCM/Test/CER-NEWS_CCCM_Aqua-FM3-MODIS-CAL-CS_RelB1_905906.20100104.hdf') #laptop
f = SD.SD('C:/Users/toha006/University/University/MSc/Models/Data/CCCM/2010/CER-NEWS_CCCM_Aqua-FM3-MODIS-CAL-CS_RelB1_905906.20100107.hdf') #Uni laptop

#view datasets
datasets_dic = f.datasets()
for idx,sds in enumerate(datasets_dic.keys()):
    print (idx,sds)
    
sds_obj = f.select('Irradiance layer center height profile') # select sds

data = sds_obj.get() # get sds data
print (data.shape) # print data dimensions

pprint.pprint( sds_obj.attributes() ) # print data attributes

"""
#read hdf5
import h5py
import os
#os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Home PC
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/') #Uni Laptop


h5f = h5py.File("2010_CCCM_tclw_av_latitude.h5", "r")

for index,key in enumerate(h5f.keys()):
    print (index, key)

"""