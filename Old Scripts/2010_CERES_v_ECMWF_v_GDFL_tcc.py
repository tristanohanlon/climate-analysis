# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 20:13:24 2019

@author: Tristan
"""

import datetime as dt  # Python standard library datetime  module
import numpy as np
from pyhdf import SD
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
from netCDF4 import MFDataset
import xarray as xr
import matplotlib.pyplot as plt

#-----------------ECMWF---------------#
# run get data xarray first

#ecmwf_lat = dataset.variables['latitude'][:] #Extract latitude data
#ecmwf_tcc = dataset.variables['tcc'][:] #Extract surface tempertaure data, keyed to time, lat and long

#ecmwf_tcloud_av = np.mean(ecmwf_tcc, axis=0) #Get mean total cloud cover over time (first key array column)
#tcloud_av1 = np.array(ecmwf_tcloud_av) # convert xarray back to numpy array

#ecmwf_lat_bins=np.arange(-89,89,1)
#ecmwf_lat_hist=np.digitize(ecmwf_lat,bins=ecmwf_lat_bins)
#ecmwf_list=[]
#for i in range(len(ecmwf_lat_bins)):
#    location=ecmwf_lat_hist==i
#    ecmwf_cf_lat=tcloud_av1[location]
#    ecmwf_nw_data=np.nanmean(ecmwf_cf_lat[(ecmwf_cf_lat<1.0)&(ecmwf_cf_lat>-0.01)])
#    ecmwf_list.append(ecmwf_nw_data)
    
#-----------------CERES---------------#
f=SD.SD('CERES/2010/CERES_2010_tcc.hdf')
#view datasets
f.datasets()

tc=f.select('tcc') 
lat=f.select('latitude (dimension)')
ceres_tc=tc.get()
lat=lat.get()

a = np.reshape(ceres_tc,180)

tcc=a/100
    
#-----------------GFDL---------------#
dataset = Dataset('GFDL/2010/clt_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc', 'r')
gfdl_lat = dataset.variables['lat'][:] #Extract latitude data
gfdl_tcc = dataset.variables['clt'][:] #Extract surface tempertaure data, keyed to time, lat and long
        
gfdl_tc = np.take(gfdl_tcc, [361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372], axis=0)
gfdl_tcloud_av = np.mean(gfdl_tc, axis=0) #Get mean surface tempertaure over time (first key array column)

gfdl_lat_bins=np.arange(-89,89,1)
gfdl_lat_hist=np.digitize(gfdl_lat,bins=gfdl_lat_bins)
gfdl_list=[]
for i in range(len(gfdl_lat_bins)):
    location=gfdl_lat_hist==i
    gfdl_cf_lat=gfdl_tcloud_av[location]
    gfdl_nw_data=np.nanmean(gfdl_cf_lat[(gfdl_cf_lat<1.0)&(gfdl_cf_lat>-0.01)])
    gfdl_list.append(gfdl_nw_data)


plt.figure()

fig, ax = plt.subplots()
ax.plot(lat,tcc, '-b', label='CERES Total Cloud Cover')
#ax.plot(ecmwf_lat_bins,ecmwf_list, '-r', label='ECMWF Total Cloud Cover')
ax.plot(gfdl_lat_bins,gfdl_list, '-g', label='GFDL Total Cloud Cover')
#ax.axis('equal')
leg = ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

plt.xlabel('Latitude')
plt.ylabel('Cloud Fraction')
plt.title('Total Cloud Fraction vs Latitude - CERES, GFDL.AM4 and ECMWF Reanalysis - 2010')
plt.grid(True)
plt.show()





















