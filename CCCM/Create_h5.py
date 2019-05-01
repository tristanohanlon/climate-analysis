# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon
"""

import h5py

# specify path and file name to create 
with h5py.File('2010_CCCM_SO_cloud_phase_temp_pressure_profile.h5', 'w') as hf:
    hf.create_dataset('Southern Ocean Cloud Fraction Profile', data=cf)
    hf.create_dataset('Southern Ocean Ice Water Content Profile', data=iw)
    hf.create_dataset('Southern Ocean Liquid Water Content Profile', data=lw)
    hf.create_dataset('Southern Ocean Temperature Profile', data=temp)
    hf.create_dataset('Southern Ocean Pressure Profile', data=pressure)
    hf.close()