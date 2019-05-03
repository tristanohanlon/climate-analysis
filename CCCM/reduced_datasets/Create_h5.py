# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon
"""

import h5py
import os

os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM')
# specify path and file name to create 
with h5py.File('2010_CCCM_SO_profile_data_test.h5', 'w') as p:
    pcf = p.create_dataset('Cloud Fraction', data=cccm_tcc_lat)
    pcf.attrs['title'] = "Cloud Fraction"
    plw = p.create_dataset('Liquid Water Content', data=cccm_tclw_lat)
    plw.attrs['title'] = "Cloud Liquid Water Content"
    plw.attrs['units'] = "gm^-3" 
    piw = p.create_dataset('Ice Water Content', data=cccm_tciw_lat)
    piw.attrs['title'] = "Cloud Ice Water Content"
    piw.attrs['units'] = "gm^-3"     
    p.close()
   
    """
import h5py

# append h5 file
f = h5py.File('2010_CCCM_SO_cloud_phase_temp_pressure_profile.h5', 'a')    
f.create.dataset('Southern Ocean Pressure Profile', data=pressure)
f.close()
"""