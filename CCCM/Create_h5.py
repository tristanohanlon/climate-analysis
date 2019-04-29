# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon
"""

import h5py

# specify path and file name to create 
with h5py.File('2010_CCCM_global_cloud_phase_lat_data2.h5', 'w') as hf:
    hf.create_dataset('tcc', data=cccm_tcc_lat)
    hf.create_dataset('tciw', data=cccm_tciw_lat)
    hf.create_dataset('tclw', data=cccm_tclw_lat)
    hf.close()