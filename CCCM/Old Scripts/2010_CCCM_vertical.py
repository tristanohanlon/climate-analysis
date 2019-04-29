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


temp=[]
alt = []
alt_c = []
counter=0

# The directory where your HDR files are stored
os.chdir('2010') 
cf = []
lw = []
iw =[]
# Load every file in the directory
for filename in os.listdir(): 
    
    # Load the file
    f = SD.SD(filename)
    alt = f.select('Irradiance layer center height profile').get().tolist() #137 levels for temp, pressure, ice and water profiles
    alt_c = f.select('Layer center height profile (clouds and aerosol)').get().tolist() #113 levels for clouds and aerosols
#    temp = temp+f.select('Temperature profile').get().tolist() #138 levels
    cf = cf+f.select('Cloud fraction profile').get().tolist() # 113 levels, #25536 values each referenced to latitude and longitude
    lw = lw+f.select('Liquid water content profile used').get().tolist() # 137 levels, #25536 values each referenced to latitude and longitude
    iw = iw+f.select('Ice water content profile used').get().tolist() # 137 levels, #25536 values each referenced to latitude and longitude


alt = np.array(alt) / 1000
cf = np.array(cf)
lw = np.array(lw)
iw = np.array(iw)

for index, item in np.ndenumerate(cf):
    if (item > 100):
        cf[index] = 0
    
for index, item in np.ndenumerate(lw):
    if (item > 20):
        lw[index] = 0
        
for index, item in np.ndenumerate(iw):
    if (item >20):
        iw[index] = 0     
        
cf = np.mean(cf, axis=0)
lw = np.mean(lw, axis=0)
iw = np.mean(iw, axis=0)


plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()

ax1.plot(cf,alt_c, '-r', label='Total Cloud Fraction')
ax2.plot(lw[25:137], alt[25:137], '-b', label='Liquid Water')
ax2.plot(iw[25:137], alt[25:137], '--b', label='Ice Water')
#ax.axis('equal')
ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);
           
ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
         ncol=4, fancybox=True, shadow=True);
           
ax1.set_xlabel('Cloud Fraction (%)', color='r')
ax1.set_ylabel('Altitude (km)')
ax2.set_xlabel('Liquid and Ice Water Content ($gm{-3}$)', color='b')

plt.title('Global Total Cloud Fraction, Ice and Water Content vs Altitude - CCCM 2010')

plt.grid(True)
plt.show()