# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon - University of Auckland & Jonathan Rogers


"""
import numpy as np
import os
from netCDF4 import Dataset
from netCDF4 import date2index
import h5py
import math
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline
from scipy import stats
import datetime
from pprint import pprint
from sklearn.impute import SimpleImputer
import constants
import time
from scipy import ndimage as nd
import matplotlib.pyplot as plt

###############################################################################
start = time.time()

# #set location
location = constants.home
# os.chdir( location + 'Data/ECMWF/' )

# #------------- Profile variables -------------#

# with Dataset( '200601-201012_ECMWF_plevel_T_cc_clw_cli_w.nc', 'r') as f:
#     p = np.flip( constants.extract_data( 'level', f ), axis = 0 )        
#     lat = np.flip( constants.extract_data( 'latitude', f ), axis = 0 )
#     ta = np.flip( np.flip( np.mean( constants.extract_data( 't', f ), axis = 0), axis = 0 ), axis = 1 ) 


# ta_p_lat = np.mean( ta, axis = -1 ) 

# #----------------------------#



# os.chdir( location + '/climate-analysis/reduced_data' )

# save_filename = 'Jan_2007_Dec_2010_ECMWF_ta_p_lat.h5'

# with h5py.File(save_filename, 'w') as h:

#     h.create_dataset('p', data=(p)) # global layer cloud liquid water fraction corresponding to ta_so
      
#     h.create_dataset('lat', data=( lat ) ) # vertical pressure velocity (Pa/s) corresponding to alt, lat

#     h.create_dataset('ta_p_lat', data= np.transpose( ta_p_lat ) ) # temperature corresponding to liq_alt and lat

#     h.close()

# end = time.time()
# print('Create dataset took:', end - start, 's')
os.chdir( location + '/climate-analysis/reduced_data' )

h5f = h5py.File( 'Jan_2007_Dec_2010_ECMWF_ta_p_lat.h5', 'r')
 
p = h5f['p'][:]
lat = h5f['lat'][:]
ta_p_lat = h5f['ta_p_lat'][:]

print(p)