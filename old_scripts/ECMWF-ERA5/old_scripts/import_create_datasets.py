# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon
"""

#---Import Latitude Variables---#

import h5py
import os
import numpy as np
os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM') #Uni Laptop
#os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CCCM') #Laptop


g = h5py.File('2010_ECMWF_global_latitude.h5', 'r')
tcc = g['tcc'][:]
tclw = g['tclw'][:]
tciw = g['tciw'][:]


b = h5py.File('2010_ECMWF_global_profile.h5', 'r')
temp = b['Temperature Profile'][:]
pressure = b['Pressure Profile'][:]
cf = b['Cloud Fraction'][:]
lw = b['Specific Liquid Water Content'][:]
iw = b['Specific Ice Water Content'][:]
cf_t = b['Cloud Fraction with Temperature'][:]
lw_t = b['Specific Liquid Water Content with Temperature'][:]
iw_t = b['Specific Ice Water Content with Temperature'][:]

c = h5py.File('2010_ECMWF_so_profile.h5', 'r')
temp_so = c['Temperature Profile'][:]
pressure_so = c['Pressure Profile'][:]
cf_so = c['Cloud Fraction'][:]
lw_so = c['Specific Liquid Water Content'][:]
iw_so = c['Specific Ice Water Content'][:]
cf_t_so = c['Cloud Fraction with Temperature'][:]
lw_t_so = c['Specific Liquid Water Content with Temperature'][:]
iw_t_so = c['Specific Ice Water Content with Temperature'][:]



###############################################################################
    

import h5py
import os
os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets') #Uni Laptop

with h5py.File('2010_ECMWF.h5', 'w') as p:
    
    p.create_dataset('tcc', data=tcc)
    p.create_dataset('tclw', data=tclw)
    p.create_dataset('tciw', data=tciw)
    
    p.create_dataset('cf', data=cf)
    p.create_dataset('cf_so', data=cf_so)
    p.create_dataset('lw', data=lw)
    p.create_dataset('lw_so', data=lw_so)
    p.create_dataset('iw', data=iw)
    p.create_dataset('iw_so', data=iw_so)

    p.create_dataset('temp', data=temp)
    p.create_dataset('temp_so', data=temp_so)
    p.create_dataset('pressure', data=pressure)
    p.create_dataset('pressure_so', data=pressure_so)

    p.create_dataset('cf_t', data=cf_t)
    p.create_dataset('cf_t_so', data=cf_t_so)
    p.create_dataset('lw_t', data=lw_t)
    p.create_dataset('lw_t_so', data=lw_t_so)
    p.create_dataset('iw_t', data=iw_t)
    p.create_dataset('iw_t_so', data=iw_t_so)

    p.close()






















