# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 14:44:25 2019

@author: Tristan O'Hanlon - University of Auckland

"""

#read hdf5
import h5py
import os
import matplotlib.pyplot as plt
import constants
import numpy as np


#--- Set Location, date period and model ---#

# specify model from the list above
model = 'CCCM' 

# specify location: home, uni, hdd, laptop
location = constants.home + 'climate-analysis/reduced_data'

os.chdir( location )


h5f = h5py.File( constants.date_cccm + '_' + model + '.h5', 'r')
cl_alt_lat = h5f['cl_alt_lat'][:]
cl_g = h5f['cl_g'][:]
cl_so = h5f['cl_so'][:]
clic_alt_lat = h5f['cli_alt_lat'][:]
cli_g = h5f['cli_g'][:]
cli_so = h5f['cli_so'][:]
clwc_alt_lat = h5f['clw_alt_lat'][:]
clw_frac_alt_lat = h5f['clw_frac_alt_lat'][:]
clw_g = h5f['clw_g'][:]
clw_so = h5f['clw_so'][:]
ta_alt_lat = h5f['ta_alt_lat'][:]
clwc_t_g = h5f['clw_t_g'][:]
clwc_t_so = h5f['clw_t_so'][:]
clw_frac_t_g = h5f['clw_frac_t_g'][:]
clw_frac_t_so = h5f['clw_frac_t_so'][:]
cl_t_g = h5f['cl_t_g'][:]
cl_t_so = h5f['cl_t_so'][:]
full_clw_alt_lat = h5f['full_clw_alt_lat'][:]
full_ta_alt_lat = h5f['full_ta_alt_lat'][:]
# clt = h5f['clt'][:]
# clt_lat_lon = h5f['clt_lat_lon'][:]
# clwvi = h5f['clwvi'][:]
# clwvi_lat_lon = h5f['clwvi_lat_lon'][:]
# clivi = h5f['clivi'][:]
# clivi_lat_lon = h5f['clivi_lat_lon'][:]
   

print(constants.globalalt_latMeanVal(cl_alt_lat[:,constants.so_idx_1:constants.so_idx_2], constants.lat[constants.so_idx_1:constants.so_idx_2]))

for index,key in enumerate(h5f.keys()):
    print (index, key)
    
