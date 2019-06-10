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
#f = SD.SD('C:/Users/toha006/University/University/MSc/Models/Data/CCCM/2006/CER-NEWS_CCCM_Aqua-FM3-MODIS-CAL-CS_RelB1_905905.20060731.hdf') #Uni laptop
f = SD.SD('E:/University/University/MSc/Models/Data/CCCM/2006/CER-NEWS_CCCM_Aqua-FM3-MODIS-CAL-CS_RelB1_905905.20060731.hdf') #Home PC


#view datasets
datasets_dic = f.datasets()
for idx,sds in enumerate(datasets_dic.keys()):
    print (idx,sds)
    
sds_obj = f.select('Ice water content profile used') # select sds

data = sds_obj.get() # get sds data
print (data.shape) # print data dimensions

pprint.pprint( sds_obj.attributes() ) # print data attributes

"""
#read hdf5
import h5py
import os
os.chdir('E:/University/University/MSc/Models/climate-analysis/gfdl/reduced_datasets') #Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/war_datasets') #Uni Laptop


h5f = h5py.File("07.2006_04.2011_gfdl.h5", "r")

for index,key in enumerate(h5f.keys()):
    print (index, key)

"""