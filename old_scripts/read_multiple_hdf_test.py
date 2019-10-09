# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: toha006
"""

import datetime as dt  # Python standard library datetime  module
import numpy as np
import os
import pylab as py
import numpy, scipy
import itertools
from scipy import interpolate
from pyhdf import SD
import matplotlib.pyplot as plt

array1=[]
array2=[]
array3=[]


count=0

os.chdir('Data/CERES/2016/') # change the working directory (where your files are stored)

for filename in os.listdir():
    try:
        print (filename)
        f = SD.SD(filename)
        subdataset_name='latitude'
        sds = f.select(subdataset_name)
        lat=sds.get()
        subdataset_name='longitude'
        sds = f.select(subdataset_name)
        lon=sds.get()
        subdataset_name='cld_amount_zon'
        sds=f.select(subdataset_name)
        data=sds.get()
        cloud_fraction=data[-1,:]
        array1=array1+cloud_fraction.ravel().tolist()
        array2=array2+lat.ravel().tolist()
        array3=array3+lon.ravel().tolist()
        count=count+1
    except:
        print (filename)
#        undefined.append(i)
                
            #instead you need to invoke a condition that simply considers the clouds that are not low clouds in a grid point and consider their effect on the average.
    cf=np.array(array1)
    lat1=np.array(array2)
    lon1=np.array(array3)
    lat_bins=np.arange(-30,5.5,0.5)[::-1]
    lon_bins=np.arange(-130,-74.5,0.5)
    lon2=[lon1[np.where(lat1==k)] for k in lat_bins]
    cf1=[cf[np.where(lat1==k)] for k in lat_bins]

n=len(os.listdir())
cf_list = list(itertools.chain.from_iterable([i]*n for i in [sum(cf[i:i+n])/n for i in range(0,len(cf),n)]))
lat_list = list(itertools.chain.from_iterable([i]*n for i in [sum(lat1[i:i+n])//n for i in range(0,len(lat1),n)]))
c = set(cf_list)
x=[float(i) for i in c]
d = set(lat1)
y=[float(i) for i in d]

plt.figure()
plt.plot(lat1 ,cf)
plt.show()
