# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

This import the reduced datasets from CCCM, ECMWF and GDFL. 
The code can select either global or southern ocean data.

"""
import time
import matplotlib.pyplot as plt
import h5py
import os

start = time.time()

# Uni Laptop
#ECMWF Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets')
2010epg = h5py.File('2010_ECMWF_global_profile_data.h5', 'r')
2010epso = h5py.File('2010_ECMWF_SO_profile_data.h5', 'r')
#CCCM Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
2010cpg = h5py.File('2010_CCCM_global_profile.h5', 'r')
2010cpso = h5py.File('2010_CCCM_SO_profile_data.h5', 'r')
2010clg = h5py.File('2010_CCCM_global_phase_lat.h5', 'r')
"""
# Home PC
#ECMWF Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets')
2010epg = h5py.File('2010_ECMWF_global_profile_data.h5', 'r')
2010epso = h5py.File('2010_ECMWF_SO_profile_data.h5', 'r')
#CCCM Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
2010cpg = h5py.File('2010_CCCM_global_profile.h5', 'r')
2010cpso = h5py.File('2010_CCCM_SO_profile_data.h5', 'r')
2010clg = h5py.File('2010_CCCM_global_phase_lat.h5', 'r')
"""





ecmwf_tcc_lat = f['Cloud Fraction Profile'][:]
ecmwf_tclw_lat = 
ecmwf_tciw_lat = 

ecmwf_tcc_alt = f['Cloud Fraction Profile'][:]
ecmwf_tclw_alt = 
ecmwf_tciw_alt = 
ecmwf_temp_alt = 
ecmwf_plevel_alt = 

ecmwf_tcc_temp = 
ecmwf_tclw_temp = 
ecmwf_tciw_temp = 

ecmwf_tcc_plevel = 
ecmwf_tclw_plevel = 
ecmwf_tciw_plevel = 

end = time.time()
print('Importing data took:', end - start, 's')
#----------------------------#

plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()
ax1.plot(ecmwf_tcc_alt[:,1],ecmwf_tcc_alt[:,0], '-r', label='Fraction Cloud Cover')
ax2.plot(ecmwf_tclw_alt[:,1],ecmwf_tclw_alt[:,0], '-b', label='Specific Cloud Liquid Water Content')
ax2.plot(ecmwf_tciw_alt[:,1],ecmwf_tciw_alt[:,0], '--b', label='Specific Cloud Ice Water Content')

#ax.axis('equal')
ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);
#ax.axis('equal')
ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
          ncol=4, fancybox=True, shadow=True);
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Specific Cloud Liquid and Ice Water Content (kg/kg) x $10^{-4}$')
ax1.set_ylabel('Altitude (km)')
plt.title(' Southern Ocean Cloud Fraction and Phase vs Altitude ECMWF 2010')
#plt.gca().invert_yaxis()

plt.grid(True)
plt.show()
