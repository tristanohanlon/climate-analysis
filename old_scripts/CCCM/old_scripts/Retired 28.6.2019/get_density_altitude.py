# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:26:33 2019

@author: Tristan O'Hanlon
"""
import h5py
import os
from scipy import integrate
import numpy as np

#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Uni Laptop
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Home PC
f = h5py.File('2010_CCCM_global_profile.h5', 'r')

temp = f['Temperature Profile'][:]
pressure = f['Pressure Profile'][:]

alt = temp[:,0]
alt = alt[1:138]
alt = alt*1000
#iw = iw[1:110]

air_density = [] #create empty list
#calculate air density at each pressure layer
air_density = (pressure[:,1] * 100) / (286.9 * temp[:,1]) * 1000
air_density = air_density[1:138]

ap = integrate.trapz(air_density, alt)
"""
lwc = lw[:,1] / air_density
iwc = iw[:,1] / air_density

temp_cf = temp[35:135]
templiwc = temp[26:135]

cf_temp = np.vstack((temp_cf[:,1],cf[:,1])).T
lwc_temp = np.vstack((templiwc[:,1],lwc)).T
iwc_temp = np.vstack((templiwc[:,1],iwc)).T

pres_cf = pressure[35:135]
presliwc = pressure[26:135]

cf_pres = np.vstack((pres_cf[:,1],cf[:,1])).T
lwc_pres = np.vstack((presliwc[:,1],lwc)).T
iwc_pres = np.vstack((presliwc[:,1],lwc)).T
"""
