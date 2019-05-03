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
#dataset = Dataset('../Data/ECMWF/pressure_levels/2010_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
dataset = Dataset('../Data/ECMWF/pressure_levels/2010_ECMWF_amon_plevels_T_cc_clw_ciw.nc', 'r')
print (dataset.variables)
tcc = dataset.variables['cc'][:]