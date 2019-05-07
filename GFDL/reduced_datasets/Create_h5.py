# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon
"""

import h5py
import os


#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GFDL') # Uni Laptop
os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL') #home PC
# specify path and file name to create 

with h5py.File('2010_gfdl_global_latitude.h5', 'w') as p:
    p.create_dataset('Cloud Fraction', data=gfdl_tcc_lat)
    p.create_dataset('Specific Liquid Water Content', data=gfdl_tclw_lat)
    p.create_dataset('Specific Ice Water Content', data=gfdl_tciw_lat)
    p.close()
    
"""
with h5py.File('2010_gfdl_global_profile.h5', 'w') as p:
    p.create_dataset('Cloud Fraction', data=cf)
    p.create_dataset('Specific Liquid Water Content', data=lw)
    p.create_dataset('Specific Ice Water Content', data=iw)
    p.create_dataset('Temperature Profile', data=temp)
    p.create_dataset('Pressure Profile', data=pressure)
    
    p.create_dataset('Cloud Fraction with Temperature', data=cf_t)
    p.create_dataset('Specific Liquid Water Content with Temperature', data=lw_t)
    p.create_dataset('Specific Ice Water Content with Temperature', data=iw_t)
    
    p.create_dataset('Cloud Fraction with Pressure', data=cf_p)
    p.create_dataset('Specific Liquid Water Content with Pressure', data=lw_p)
    p.create_dataset('Specific Ice Water Content with Pressure', data=iw_p)
    p.close()
"""

"""
with h5py.File('2010_gfdl_so_profile.h5', 'w') as p:
    p.create_dataset('Cloud Fraction', data=cf)
    p.create_dataset('Specific Liquid Water Content', data=lw)
    p.create_dataset('Specific Ice Water Content', data=iw)
    p.create_dataset('Temperature Profile', data=temp)
    p.create_dataset('Pressure Profile', data=pressure)
    
    p.create_dataset('Cloud Fraction with Temperature', data=cf_t)
    p.create_dataset('Specific Liquid Water Content with Temperature', data=lw_t)
    p.create_dataset('Specific Ice Water Content with Temperature', data=iw_t)
    
    p.create_dataset('Cloud Fraction with Pressure', data=cf_p)
    p.create_dataset('Specific Liquid Water Content with Pressure', data=lw_p)
    p.create_dataset('Specific Ice Water Content with Pressure', data=iw_p)
    p.close()
"""      