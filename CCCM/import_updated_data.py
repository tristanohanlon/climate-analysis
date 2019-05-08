# -*- coding: utf-8 -*-
"""
Created on Wed May  8 11:56:14 2019

@author: Tristan O'Hanlon

Import updated CCCM data 
"""
import matplotlib.pyplot as plt
import h5py
import os

#os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Home PC
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Uni Laptop

f = h5py.File('2010_CCCM_global_latitude.h5', 'r')
cf = f['Cloud Fraction'][:]

#os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Home PC
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/') #Uni Laptop

p = h5py.File('2010_CCCM_tciw_lat.h5', 'r')
iw = p['Specific Ice Water Content'][:]

#os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets') #Home PC
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/') #Uni Laptop

p = h5py.File('2010_CCCM_tclw_lat.h5', 'r')
lw = p['Specific Liquid Water Content'][:]
