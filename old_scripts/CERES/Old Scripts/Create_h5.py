# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon
"""

import h5py

# specify path and file name to create 
with h5py.File('2010_CERES_global_cloud_phase_profile_data.h5', 'w') as hf:
    hf.create_dataset('Cloud Fraction Profile', data=ceres_tcc_alt)
    hf.create_dataset('Liquid Water Content Profile', data=ceres_tclw_alt)
    hf.create_dataset('Ice Water Content Profile', data=ceres_tciw_alt)
    hf.close()    

with h5py.File('2010_CERES_global_cloud_phase_lat_data.h5', 'w') as hf:
    hf.create_dataset('Total Cloud Cover', data=ceres_tcc_lat)
    hf.create_dataset('Total Liquid Water Cloud Cover', data=tclw)
    hf.create_dataset('Total Ice Water Cloud Cover', data=tciw)
    hf.close()  