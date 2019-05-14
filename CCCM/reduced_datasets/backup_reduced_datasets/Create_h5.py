# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon
"""

import h5py
import os
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM') # Home PC
#os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM') #Uni Laptop

# specify path and file name to create 
with h5py.File('2010_CCCM_global_latitude.h5', 'w') as p:
    p.create_dataset('Cloud Fraction', data=cccm_tcc_lat)
    p.create_dataset('Liquid Water Content', data=cccm_tclw_lat)
    p.create_dataset('Ice Water Content', data=cccm_tciw_lat)
#    p.create_dataset('Mean Specific Liquid Water Content', data=cccm_tclw_av_lat)
#    p.create_dataset('Mean Specific Ice Water Content', data=cccm_tciw_av_lat)
    p.close()
""" 
import h5py
import os
#os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM') # Home PC
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM') #Uni Laptop

with h5py.File('2010_CCCM_global_profile.h5', 'w') as p:
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
import h5py
import os
#os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM') # Home PC
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM') #Uni Laptop

with h5py.File('2010_CCCM_so_profile.h5', 'w') as p:
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
import h5py
import os
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM') # Home PC
# specify path and file name to create 
with h5py.File('2010_lat_raw.h5', 'w') as p:
    p.create_dataset('Latitude', data=lat)
    p.create_dataset('Temperature Profile', data=temp)
    p.create_dataset('Temperature Profile', data=pressure)
    p.create_dataset('Liquid Water Content', data=tclw)
    p.create_dataset('Ice Water Content', data=tciw)
    p.close()
""" 