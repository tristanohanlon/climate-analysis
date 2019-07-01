# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon
"""

import datetime as dt  # Python standard library datetime  module
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Specify nc file
#dataset = Dataset('E:/University/University/MSc/Models/Data/CMIP6/cesm2.1_cam6/f.e21.FHIST_BGC.f09_f09_mg17.CMIP6-AMIP.001.cam.h0.CLOUD.195001-201412.nc', 'r') #Home PC
#dataset = Dataset('D:/MSc/Models/Data/CALIPSO-GOCCP/3D_CloudFraction_Phase330m_200606-201803_avg_CFMIP2_sat_3.1.2.nc', 'r') #ext HDD
dataset = Dataset('//Synthesis/E/University/University/MSc/Models/Data/CMIP6/gfdl_am4/clw_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r') #Laptop
#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/CMIP5/gfdl_cm3_rcp4.5/cl_Amon_GFDL-CM3_rcp45_r1i1p1_200601-201012.nc', 'r') #Uni Laptop

print (dataset.variables)
#x = dataset.variables['ta'][:]