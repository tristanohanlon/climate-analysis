# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import datetime as dt  # Python standard library datetime  module
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
from netCDF4 import MFDataset
import xarray as xr
import matplotlib.pyplot as plt

#dataset = Dataset('2010/2010_tcc.nc', 'r')
dataset = xr.open_mfdataset('2010/2010_ECMWF_tcc.nc', chunks={'time': 200})
#elon = dataset.variables['longitude'][:] #Extract longitude data
elat = dataset.variables['latitude'][:] #Extract latitude data
etcc = dataset.variables['tcc'][:] #Extract total cloud cover, keyed to time, lat and long