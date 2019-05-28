# -*- coding: utf-8 -*-
"""
Created on Mon May  6 08:37:25 2019

@author: Tristan O'Hanlon

Convert total column cloud liquid water (LWP) (kgm^-2) and total column cloud ice water (IWP) (kgm^-2)
to specific ice and water content (kg/kg)

Find the total colomn air (AP) (kgm^-2) by summing air density at all pressure levels in profile data
then divideLWP and IWP by AP.

From ideal gas law:
density = pressure / (R * temperature)
R = 286.9 J/Kg.K
"""

import matplotlib.pyplot as plt
import h5py
import os
import numpy as np
from scipy import integrate

#os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Home PC
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets') #Uni Laptop

p = h5py.File('2010_ECMWF_global_latitude.h5', 'r')
f = h5py.File('2010_ECMWF_global_profile.h5', 'r')


lw = f['Specific Liquid Water Content'][:]
iw = f['Specific Ice Water Content'][:]
temp = f['Temperature Profile'][:]
pressure = f['Pressure Profile'][:]

ecmwf_tcc_lat = p['Cloud Fraction'][:]
ecmwf_tclw_lat = p['Specific Ice Water Content'][:]
ecmwf_tciw_lat = p['Specific Ice Water Content'][:]

air_density = [] #create empty list
#calculate air density at each pressure layer
air_density = (pressure[:,1] * 100) / (286.9 * temp[:,1])

#mean air density:
air_density_mean = np.mean(air_density)

#mean specific cloud liquid and ice water density
lw_av = ecmwf_tclw_lat[:,1] / (np.amax(lw[:,0])*1000)
iw_av = ecmwf_tciw_lat[:,1] / (np.amax(iw[:,0])*1000)

#mean specific cloud liquid and ice water content
ecmwf_tclw_av_lat = lw_av / air_density_mean
ecmwf_tciw_av_lat = iw_av / air_density_mean

#integrate air density with altitude to get total air column
#ap = integrate.trapz(air_density, ecmwf_plevel_alt_g[:,0]*1000)

#convert LWP to specific liquid water content
#ecmwf_tclw_lat = ecmwf_tclw_lat / -ap

#convert IWP to specific ice water content
#ecmwf_tciw_lat = ecmwf_tciw_lat / -ap

