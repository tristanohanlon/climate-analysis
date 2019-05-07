# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon
"""

import h5py
import os

# Uni Laptop
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/ECMWF')
#home PC
os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF')
# specify path and file name to create 
"""
with h5py.File('2010_ECMWF_so_profile.h5', 'w') as p:
    p.create_dataset('Cloud Fraction', data=ecmwf_tcc_alt)
    p.create_dataset('Specific Liquid Water Content', data=ecmwf_tclw_alt)
    p.create_dataset('Specific Ice Water Content', data=ecmwf_tciw_alt)
    p.create_dataset('Temperature Profile', data=ecmwf_temp_alt)
    p.create_dataset('Pressure Profile', data=ecmwf_plevel_alt)
    
    p.create_dataset('Cloud Fraction with Temperature', data=ecmwf_tcc_temp)
    p.create_dataset('Specific Liquid Water Content with Temperature', data=ecmwf_tclw_temp)
    p.create_dataset('Specific Ice Water Content with Temperature', data=ecmwf_tciw_temp)
    
    p.create_dataset('Cloud Fraction with Pressure', data=ecmwf_tcc_plevel)
    p.create_dataset('Specific Liquid Water Content with Pressure', data=ecmwf_tclw_plevel)
    p.create_dataset('Specific Ice Water Content with Pressure', data=ecmwf_tciw_plevel)
    p.close()
"""   


with h5py.File('2010_ECMWF_global_latitude.h5', 'w') as p:
    p.create_dataset('Cloud Fraction', data=ecmwf_tcc_lat)
    p.create_dataset('Specific Liquid Water Content', data=ecmwf_tclw_lat)
    p.create_dataset('Specific Ice Water Content', data=ecmwf_tciw_lat)

    p.close()