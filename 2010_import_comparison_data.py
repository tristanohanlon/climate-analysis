# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

This import the reduced datasets from 2010 CCCM, ECMWF and GDFL. 
The code can select either global or southern ocean data.

"""
import time
import matplotlib.pyplot as plt
import h5py
import os

start = time.time()


#---Importing Data from Reduced Datasets---#
"""
# Uni Laptop
#ECMWF Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets')
epg = h5py.File('2010_ECMWF_global_profile.h5', 'r')
epso = h5py.File('2010_ECMWF_SO_profile.h5', 'r')
elg = h5py.File('2010_ECMWF_global_latitude.h5', 'r')
#CCCM Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
cpg = h5py.File('2010_CCCM_global_profile.h5', 'r')
cpso = h5py.File('2010_CCCM_SO_profile.h5', 'r')
clg = h5py.File('2010_CCCM_global_latitude.h5', 'r')
#GDFL Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GDFL/reduced_datasets')
gpg = h5py.File('2010_GDFL_global_profile.h5', 'r')
gpso = h5py.File('2010_GDFL_SO_profile.h5', 'r')
glg = h5py.File('2010_GDFL_global_latitude.h5', 'r')
"""
# Home PC
#ECMWF Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets')
epg = h5py.File('2010_ECMWF_global_profile.h5', 'r')
epso = h5py.File('2010_ECMWF_SO_profile.h5', 'r')
elg = h5py.File('2010_ECMWF_global_latitude.h5', 'r')

#CCCM Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
cpg = h5py.File('2010_CCCM_global_profile.h5', 'r')
cpso = h5py.File('2010_CCCM_SO_profile.h5', 'r')
clg = h5py.File('2010_CCCM_global_latitude.h5', 'r')
#GDFL Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/GDFL/reduced_datasets')
gpg = h5py.File('2010_GDFL_global_profile.h5', 'r')
gpso = h5py.File('2010_GDFL_SO_profile.h5', 'r')
glg = h5py.File('2010_GDFL_global_latitude.h5', 'r')

############################################################################### ECMWF Data

#---ECMWF Global Latitude Data---#

ecmwf_tcc_lat_g = elg['Cloud Fraction'][:] # 0-1
ecmwf_tclw_lat_g = elg['Specific Liquid Water Content'][:] #kg/kg
ecmwf_tciw_lat_g = elg['Specific Ice Water Content'][:] #kg/kg

#---ECMWF Southern Ocean Latitude Data---#

ecmwf_tcc_lat_so = ecmwf_tcc_lat_g # 0-1
ecmwf_tclw_lat_so = ecmwf_tclw_lat_g #kg/kg
ecmwf_tciw_lat_so = ecmwf_tciw_lat_g #kg/kg

#---ECMWF Global Profile---#

ecmwf_tcc_alt_g = epg['Cloud Fraction'][:] # 0-1
ecmwf_tclw_alt_g = epg['Specific Liquid Water Content'][:] #kg/kg
ecmwf_tciw_alt_g = epg['Specific Ice Water Content'][:] #kg/kg
ecmwf_temp_alt_g = epg['Temperature Profile'][:] #K
ecmwf_plevel_alt_g = epg['Pressure Profile'][:] #hPa

ecmwf_tcc_temp_g = epg['Cloud Fraction with Temperature'][:] # 0-1
ecmwf_tclw_temp_g = epg['Specific Liquid Water Content with Temperature'][:] #kg/kg
ecmwf_tciw_temp_g = epg['Specific Ice Water Content with Temperature'][:] #kg/kg

ecmwf_tcc_plevel_g = epg['Cloud Fraction with Pressure'][:] # 0-1
ecmwf_tclw_plevel_g = epg['Specific Liquid Water Content with Pressure'][:] #kg/kg
ecmwf_tciw_plevel_g = epg['Specific Ice Water Content with Pressure'][:] #kg/kg

#---ECMWF Southern Ocean Profile---#

ecmwf_tcc_alt_so = epso['Cloud Fraction'][:] # 0-1
ecmwf_tclw_alt_so = epso['Specific Liquid Water Content'][:] #kg/kg
ecmwf_tciw_alt_so = epso['Specific Ice Water Content'][:] #kg/kg
ecmwf_temp_alt_so = epso['Temperature Profile'][:] #K
ecmwf_plevel_alt_so = epso['Pressure Profile'][:] #hPa

ecmwf_tcc_temp_so = epso['Cloud Fraction with Temperature'][:] # 0-1
ecmwf_tclw_temp_so = epso['Specific Liquid Water Content with Temperature'][:] #kg/kg
ecmwf_tciw_temp_so = epso['Specific Ice Water Content with Temperature'][:] #kg/kg

ecmwf_tcc_plevel_so = epso['Cloud Fraction with Pressure'][:] # 0-1
ecmwf_tclw_plevel_so = epso['Specific Liquid Water Content with Pressure'][:] #kg/kg
ecmwf_tciw_plevel_so = epso['Specific Ice Water Content with Pressure'][:] #kg/kg

############################################################################### CCCM Data

#---CCCM Global Latitude Data---#

cccm_tcc_lat_g = clg['Cloud Fraction'][:] # 0-1
cccm_tclw_lat_g = clg['Specific Liquid Water Content'][:] #kg/kg
cccm_tciw_lat_g = clg['Specific Ice Water Content'][:] #kg/kg

#---CCCM Southern Ocean Latitude Data---#

cccm_tcc_lat_so = cccm_tcc_lat_g # 0-1
cccm_tclw_lat_so = cccm_tclw_lat_g #kg/kg
cccm_tciw_lat_so = cccm_tciw_lat_g #kg/kg

#---CCCM Global Profile---#

cccm_tcc_alt_g = cpg['Cloud Fraction'][:] # 0-1
cccm_tclw_alt_g = cpg['Specific Liquid Water Content'][:] #kg/kg
cccm_tciw_alt_g = cpg['Specific Ice Water Content'][:] #kg/kg
cccm_temp_alt_g = cpg['Temperature Profile'][:] #K
cccm_plevel_alt_g = cpg['Pressure Profile'][:] #hPa

cccm_tcc_temp_g = cpg['Cloud Fraction with Temperature'][:] # 0-1
cccm_tclw_temp_g = cpg['Specific Liquid Water Content with Temperature'][:] #kg/kg
cccm_tciw_temp_g = cpg['Specific Ice Water Content with Temperature'][:] #kg/kg

cccm_tcc_plevel_g = cpg['Cloud Fraction with Pressure'][:] # 0-1
cccm_tclw_plevel_g = cpg['Specific Liquid Water Content with Pressure'][:] #kg/kg
cccm_tciw_plevel_g = cpg['Specific Ice Water Content with Pressure'][:] #kg/kg

#---CCCM Southern Ocean Profile---#

cccm_tcc_alt_so = cpso['Cloud Fraction'][:] # 0-1
cccm_tclw_alt_so = cpso['Specific Liquid Water Content'][:] #kg/kg
cccm_tciw_alt_so = cpso['Specific Ice Water Content'][:] #kg/kg
cccm_temp_alt_so = cpso['Temperature Profile'][:] #K
cccm_plevel_alt_so = cpso['Pressure Profile'][:] #hPa

cccm_tcc_temp_so = cpso['Cloud Fraction with Temperature'][:] # 0-1
cccm_tclw_temp_so = cpso['Specific Liquid Water Content with Temperature'][:] #kg/kg
cccm_tciw_temp_so = cpso['Specific Ice Water Content with Temperature'][:] #kg/kg

cccm_tcc_plevel_so = cpso['Cloud Fraction with Pressure'][:] # 0-1
cccm_tclw_plevel_so = cpso['Specific Liquid Water Content with Pressure'][:] #kg/kg
cccm_tciw_plevel_so = cpso['Specific Ice Water Content with Pressure'][:] #kg/kg

############################################################################### GDFL Data

#---GDFL Global Latitude Data---#

gdfl_tcc_lat_g = glg['Cloud Fraction'][:] # 0-1
gdfl_tclw_lat_g = glg['Specific Liquid Water Content'][:] #kg/kg
gdfl_tciw_lat_g = glg['Specific Ice Water Content'][:] #kg/kg

#---GDFL Southern Ocean Latitude Data---#

gdfl_tcc_lat_so = gdfl_tcc_lat_g # 0-1
gdfl_tclw_lat_so = gdfl_tclw_lat_g #kg/kg
gdfl_tciw_lat_so = gdfl_tciw_lat_g #kg/kg

#---GDFL Global Profile---#

gdfl_tcc_alt_g = gpg['Cloud Fraction'][:] # 0-1
gdfl_tclw_alt_g = gpg['Specific Liquid Water Content'][:] #kg/kg
gdfl_tciw_alt_g = gpg['Specific Ice Water Content'][:] #kg/kg
gdfl_temp_alt_g = gpg['Temperature Profile'][:] #K
gdfl_plevel_alt_g = gpg['Pressure Profile'][:] #hPa

gdfl_tcc_temp_g = gpg['Cloud Fraction with Temperature'][:] # 0-1
gdfl_tclw_temp_g = gpg['Specific Liquid Water Content with Temperature'][:] #kg/kg
gdfl_tciw_temp_g = gpg['Specific Ice Water Content with Temperature'][:] #kg/kg

gdfl_tcc_plevel_g = gpg['Cloud Fraction with Pressure'][:] # 0-1
gdfl_tclw_plevel_g = gpg['Specific Liquid Water Content with Pressure'][:] #kg/kg
gdfl_tciw_plevel_g = gpg['Specific Ice Water Content with Pressure'][:] #kg/kg

#---GDFL Southern Ocean Profile---#

gdfl_tcc_alt_so = gpso['Cloud Fraction'][:] # 0-1
gdfl_tclw_alt_so = gpso['Specific Liquid Water Content'][:] #kg/kg
gdfl_tciw_alt_so = gpso['Specific Ice Water Content'][:] #kg/kg
gdfl_temp_alt_so = gpso['Temperature Profile'][:] #K
gdfl_plevel_alt_so = gpso['Pressure Profile'][:] #hPa

gdfl_tcc_temp_so = gpso['Cloud Fraction with Temperature'][:] # 0-1
gdfl_tclw_temp_so = gpso['Specific Liquid Water Content with Temperature'][:] #kg/kg
gdfl_tciw_temp_so = gpso['Specific Ice Water Content with Temperature'][:] #kg/kg

gdfl_tcc_plevel_so = gpso['Cloud Fraction with Pressure'][:] # 0-1
gdfl_tclw_plevel_so = gpso['Specific Liquid Water Content with Pressure'][:] #kg/kg
gdfl_tciw_plevel_so = gpso['Specific Ice Water Content with Pressure'][:] #kg/kg



############################################################################### End importing Data

end = time.time()
print('Importing data took:', end - start, 's')

############################################################################### Global Latitude Plots

#---Plot Global Cloud Fraction with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tcc_lat_g[:,0],cccm_tcc_lat_g[:,1], '-r', label='CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tcc_lat_g[:,0],ecmwf_tcc_lat_g[:,1], '-b', label='ECMWF.ERA5 Reanalysis Model')
ax.plot(gdfl_tcc_lat_g[:,0],gdfl_tcc_lat_g[:,1], '-g', label='GDFL.AM4 Model')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Fraction')

plt.title('2010 Global Cloud Fraction vs Latitude')

plt.grid(True)
plt.show()
"""

#---Plot Global Specific Liquid Water Content with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_lat_g[:,0],cccm_tclw_lat_g[:,1], '-r', label='CCCM - Satellite')
ax.plot(ecmwf_tclw_lat_g[:,0],ecmwf_tclw_lat_g[:,1], '-b', label='ECMWF.ERA5 Reanalysis - Model')
ax.plot(gdfl_tclw_lat_g[:,0],gdfl_tclw_lat_g[:,1], '-g', label='GDFL.AM4 - Model')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Latitude')
ax.set_ylabel('Specific Liquid Water Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Global Specific Liquid Water Content vs Latitude')

plt.grid(True)
plt.show()
"""

#---Plot Global Specific Ice Water Content with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tciw_lat_g[:,0],cccm_tciw_lat_g[:,1], '--r', label='CCCM - Satellite')
ax.plot(ecmwf_tciw_lat_g[:,0],ecmwf_tciw_lat_g[:,1], '--b', label='ECMWF.ERA5 Reanalysis - Model')
ax.plot(gdfl_tciw_lat_g[:,0],gdfl_tciw_lat_g[:,1], '--g', label='GDFL.AM4 - Model')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Latitude')
ax.set_ylabel('Specific Ice Water Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Global Specific Ice Water Content vs Latitude')

plt.grid(True)
plt.show()
"""

#---Plot Global Comparison Specific Liquid and Ice Water Content with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_lat_g[:,0],cccm_tclw_lat_g[:,1], '-r', label='Liquid - CCCM - Satellite')
ax.plot(ecmwf_tclw_lat_g[:,0],ecmwf_tclw_lat_g[:,1], '-b', label='Liquid - ECMWF.ERA5 Reanalysis - Model')
ax.plot(gdfl_tclw_lat_g[:,0],gdfl_tclw_lat_g[:,1], '-g', label='Liquid - GDFL.AM4 - Model')
ax.plot(cccm_tciw_lat_g[:,0],cccm_tciw_lat_g[:,1], '--r', label='Ice - CCCM - Satellite')
ax.plot(ecmwf_tciw_lat_g[:,0],ecmwf_tciw_lat_g[:,1], '--b', label='Ice - ECMWF.ERA5 Reanalysis - Model')
ax.plot(gdfl_tciw_lat_g[:,0],gdfl_tciw_lat_g[:,1], '--g', label='Ice - GDFL.AM4 - Model')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Latitude')
ax.set_ylabel('Specific Liquid and Ice Water Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Global Specific Liquid and Ice Water Content vs Latitude')

plt.grid(True)
plt.show()
"""

############################################################################### Southern Ocean Latitude Plots

#---Plot Southern Ocean Cloud Fraction with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tcc_lat_so[:,0],cccm_tcc_lat_so[:,1], '-r', label='CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tcc_lat_so[:,0],ecmwf_tcc_lat_so[:,1], '-b', label='ECMWF.ERA5 Reanalysis Model')
ax.plot(gdfl_tcc_lat_so[:,0],gdfl_tcc_lat_so[:,1], '-g', label='GDFL.AM4 Model')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Fraction')

plt.title('2010 Southern Ocean Cloud Fraction vs Latitude')

plt.grid(True)
plt.show()
"""

#---Plot Southern Ocean Specific Liquid Water Content with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_lat_so[:,0],cccm_tclw_lat_so[:,1], '-r', label='CCCM - Satellite')
ax.plot(ecmwf_tclw_lat_so[:,0],ecmwf_tclw_lat_so[:,1], '-b', label='ECMWF.ERA5 Reanalysis - Model')
ax.plot(gdfl_tclw_lat_so[:,0],gdfl_tclw_lat_so[:,1], '-g', label='GDFL.AM4 - Model')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Latitude')
ax.set_ylabel('Specific Liquid Water Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Southern Ocean Specific Liquid Water Content vs Latitude')

plt.grid(True)
plt.show()
"""

#---Plot Southern Ocean Specific Ice Water Content with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tciw_lat_so[:,0],cccm_tciw_lat_so[:,1], '--r', label='CCCM - Satellite')
ax.plot(ecmwf_tciw_lat_so[:,0],ecmwf_tciw_lat_so[:,1], '--b', label='ECMWF.ERA5 Reanalysis - Model')
ax.plot(gdfl_tciw_lat_so[:,0],gdfl_tciw_lat_so[:,1], '--g', label='GDFL.AM4 - Model')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Latitude')
ax.set_ylabel('Specific Ice Water Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Southern Ocean Specific Ice Water Content vs Latitude')

plt.grid(True)
plt.show()
"""

#---Plot Southern Ocean Comparison Specific Liquid and Ice Water Content with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_lat_so[:,0],cccm_tclw_lat_so[:,1], '-r', label='Liquid - CCCM - Satellite')
ax.plot(ecmwf_tclw_lat_so[:,0],ecmwf_tclw_lat_so[:,1], '-b', label='Liquid - ECMWF.ERA5 Reanalysis - Model')
ax.plot(gdfl_tclw_lat_so[:,0],gdfl_tclw_lat_so[:,1], '-g', label='Liquid - GDFL.AM4 - Model')
ax.plot(cccm_tciw_lat_so[:,0],cccm_tciw_lat_so[:,1], '--r', label='Ice - CCCM - Satellite')
ax.plot(ecmwf_tciw_lat_so[:,0],ecmwf_tciw_lat_so[:,1], '--b', label='Ice - ECMWF.ERA5 Reanalysis - Model')
ax.plot(gdfl_tciw_lat_so[:,0],gdfl_tciw_lat_so[:,1], '--g', label='Ice - GDFL.AM4 - Model')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Latitude')
ax.set_ylabel('Specific Liquid and Ice Water Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Southern Ocean Specific Liquid and Ice Water Content vs Latitude')

plt.grid(True)
plt.show()
"""




############################################################################### Global Altitude Plots

#---Plot Global Cloud Fraction Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tcc_alt_g[:,0],cccm_tcc_alt_g[:,1], '-r', label='CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tcc_alt_g[:,0],ecmwf_tcc_alt_g[:,1], '-b', label='ECMWF.ERA5 Reanalysis Model')
ax.plot(gdfl_tcc_alt_g[:,0],gdfl_tcc_alt_g[:,1], '-g', label='GDFL.AM4 Model')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Altitude (km)')
ax.set_ylabel('Cloud Fraction')

plt.title('2010 Global Cloud Fraction vs Altitude')

plt.grid(True)
plt.show()
"""

#---Plot Global Specific Liquid Water Content Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_alt_g[:,0],cccm_tclw_alt_g[:,1], '-r', label='CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tclw_alt_g[:,0],ecmwf_tclw_alt_g[:,1], '-b', label='ECMWF.ERA5 Reanalysis Model')
ax.plot(gdfl_tclw_alt_g[:,0],gdfl_tclw_alt_g[:,1], '-g', label='GDFL.AM4 Model')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Altitude (km)')
ax.set_ylabel('Specific Liquid Water Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Global Specific Liquid Water Content vs Altitude')

plt.grid(True)
plt.show()
"""

#---Plot Global Specific Ice Water Content Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tciw_alt_g[:,0],cccm_tciw_alt_g[:,1], '--r', label='CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tciw_alt_g[:,0],ecmwf_tciw_alt_g[:,1], '--b', label='ECMWF.ERA5 Reanalysis Model')
ax.plot(gdfl_tciw_alt_g[:,0],gdfl_tciw_alt_g[:,1], '--g', label='GDFL.AM4 Model')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Altitude (km)')
ax.set_ylabel('Specific Ice Water Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Global Specific Ice Water Content vs Altitude')

plt.grid(True)
plt.show()
"""

#---Plot Global Comparison Specific Liquid and Ice Water Content Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_alt_g[:,0],cccm_tclw_alt_g[:,1], '-r', label='Liquid - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tclw_alt_g[:,0],ecmwf_tclw_alt_g[:,1], '-b', label='Liquid - ECMWF.ERA5 Reanalysis Model')
ax.plot(gdfl_tclw_alt_g[:,0],gdfl_tclw_alt_g[:,1], '-g', label='Liquid - GDFL.AM4 Model')
ax.plot(cccm_tciw_alt_g[:,0],cccm_tciw_alt_g[:,1], '--r', label='Ice - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tciw_alt_g[:,0],ecmwf_tciw_alt_g[:,1], '--b', label='Ice - ECMWF.ERA5 Reanalysis Model')
ax.plot(gdfl_tciw_alt_g[:,0],gdfl_tciw_alt_g[:,1], '--g', label='Ice - GDFL.AM4 Model')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Altitude (km)')
ax.set_ylabel('Specific Liquid and Ice Water Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Global Specific Liquid and Ice Water Content vs Altitude')

plt.grid(True)
plt.show()
"""

#---Plot Global Temperature Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_temp_alt_g[:,0],cccm_temp_alt_g[:,1], '-r', label='CCCM Merged Satellite Dataset')
ax.plot(ecmwf_temp_alt_g[:,0],ecmwf_temp_alt_g[:,1], '-b', label='ECMWF.ERA5 Reanalysis Model')
ax.plot(gdfl_temp_alt_g[:,0],gdfl_temp_alt_g[:,1], '-g', label='GDFL.AM4 Model')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Altitude (km)')
ax.set_ylabel('Temperature (K)')

plt.title('2010 Global Temperature vs Altitude')

plt.grid(True)
plt.show()
"""

#---Plot Global Pressure Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_plevel_alt_g[:,0],cccm_plevel_alt_g[:,1], '-r', label='CCCM Merged Satellite Dataset')
ax.plot(ecmwf_plevel_alt_g[:,0],ecmwf_plevel_alt_g[:,1], '-b', label='ECMWF.ERA5 Reanalysis Model')
ax.plot(gdfl_plevel_alt_g[:,0],gdfl_plevel_alt_g[:,1], '-g', label='GDFL.AM4 Model')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Altitude (km)')
ax.set_ylabel('Pressure (hPa)')

plt.title('2010 Global Temperature vs Altitude')
plt.gca().invert_yaxis()
plt.grid(True)
plt.show()
"""

#---Plot Global Specific Liquid Water Content with Temperature---#

#---Plot Global Specific Ice Water Content with Presure---#

#---Plot Global Presure with Temperature---#

############################################################################### Southern Ocean Altitude Plots

#---Plot Southern Ocean Cloud Fraction Altitude Profile---#

#---Plot Southern Ocean Specific Liquid Water Content Altitude Profile---#

#---Plot Southern Ocean Specific Ice Water Content Altitude Profile---#

#---Plot Southern Ocean Temperature Altitude Profile---#

#---Plot Southern Ocean Pressure Altitude Profile---#

#---Plot Southern Ocean Specific Liquid Water Content with Temperature---#

#---Plot Southern Ocean Specific Ice Water Content with Presure---#

#---Plot Southern Ocean Presure with Temperature---#








