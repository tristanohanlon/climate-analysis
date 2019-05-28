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
dataset = Dataset('E:/University/University/MSc/Models/Data/CAM5/tclw/b.e11.BRCP45C5CNBDRD.f09_g16.001.cam.h0.CLDLIQ.200601-208012.nc', 'r') #ext HDD
#dataset = Dataset('D:/MSc/Models/Data/CESM1(CAM5)/ps/b.e11.BRCP45C5CNBDRD.f09_g16.001.cam.h0.PS.200601-208012.nc', 'r') #ext HDD
#dataset = Dataset('C:/Users/tristan/University/University/MSc/Models/Data/GFDL/clt_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r') #Laptop
#dataset = Dataset('C:/Users/toha006/University/University/MSc/Models/Data/GFDL/cli_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r') #Uni Laptop

print (dataset.variables)
#tcc = dataset.variables['tclw'][:]