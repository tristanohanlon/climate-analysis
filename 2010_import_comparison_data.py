# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

This import the reduced datasets from 2010 CCCM, ECMWF and GFDL. 
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
epso = h5py.File('2010_ECMWF_so_profile.h5', 'r')
elg = h5py.File('2010_ECMWF_global_latitude.h5', 'r')

#CCCM Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
c = h5py.File('2010_CCCM.h5', 'r')

#GFDL Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/gfdl/reduced_datasets')
gpg = h5py.File('2010_gfdl_global_profile.h5', 'r')
gpso = h5py.File('2010_gfdl_so_profile.h5', 'r')
glg = h5py.File('2010_gfdl_global_latitude.h5', 'r')
"""

# Home PC
#ECMWF Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets')
epg = h5py.File('2010_ECMWF_global_profile.h5', 'r')
epso = h5py.File('2010_ECMWF_so_profile.h5', 'r')
elg = h5py.File('2010_ECMWF_global_latitude.h5', 'r')

#CCCM Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
c = h5py.File('2010_CCCM.h5', 'r')

#GFDL Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL/reduced_datasets')
gpg = h5py.File('2010_gfdl_global_profile.h5', 'r')
gpso = h5py.File('2010_gfdl_so_profile.h5', 'r')
glg = h5py.File('2010_gfdl_global_latitude.h5', 'r')

"""
# Laptop
#ECMWF Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets')
epg = h5py.File('2010_ECMWF_global_profile.h5', 'r')
epso = h5py.File('2010_ECMWF_so_profile.h5', 'r')
#elg = h5py.File('2010_ECMWF_global_latitude.h5', 'r')

#CCCM Data
#os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
c = h5py.File('2010_CCCM.h5', 'r')

#GFDL Data
#os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/gfdl/reduced_datasets')
gpg = h5py.File('2010_gfdl_global_profile.h5', 'r')
gpso = h5py.File('2010_gfdl_so_profile.h5', 'r')
glg = h5py.File('2010_gfdl_global_latitude.h5', 'r')
"""

############################################################################### ECMWF Data

#---ECMWF Global Latitude Data---#

ecmwf_tcc_lat_g = elg['Cloud Fraction'][:] # 0-1
ecmwf_tclw_lat_g = elg['Liquid Water Content'][:] #kg/kg
ecmwf_tciw_lat_g = elg['Ice Water Content'][:] #kg/kg

#---ECMWF Southern Ocean Latitude Data---#

ecmwf_tcc_lat_so = ecmwf_tcc_lat_g[ecmwf_tcc_lat_g[:,0]>=-70]
ecmwf_tcc_lat_so = ecmwf_tcc_lat_so[ecmwf_tcc_lat_so[:,0]<=-50] # 0-1

ecmwf_tclw_lat_so = ecmwf_tclw_lat_g[ecmwf_tclw_lat_g[:,0]>=-70]
ecmwf_tclw_lat_so = ecmwf_tclw_lat_so[ecmwf_tclw_lat_so[:,0]<=-50] #kg/kg

ecmwf_tciw_lat_so = ecmwf_tciw_lat_g[ecmwf_tciw_lat_g[:,0]>=-70]
ecmwf_tciw_lat_so = ecmwf_tciw_lat_so[ecmwf_tciw_lat_so[:,0]<=-50] #kg/kg

#---ECMWF Global Profile---#

ecmwf_tcc_alt_g = epg['Cloud Fraction'][8:37] # 0-1
ecmwf_tclw_alt_g = epg['Specific Liquid Water Content'][14:37] #kg/kg
ecmwf_tciw_alt_g = epg['Specific Ice Water Content'][9:37] #kg/kg
ecmwf_temp_alt_g = epg['Temperature Profile'][:] #K
ecmwf_plevel_alt_g = epg['Pressure Profile'][:] #hPa

ecmwf_tcc_temp_g = epg['Cloud Fraction with Temperature'][13:37] # 0-1
ecmwf_tclw_temp_g = epg['Specific Liquid Water Content with Temperature'][17:37] #kg/kg
ecmwf_tciw_temp_g = epg['Specific Ice Water Content with Temperature'][10:37] #kg/kg

ecmwf_tcc_plevel_g = epg['Cloud Fraction with Pressure'][8:37] # 0-1
ecmwf_tclw_plevel_g = epg['Specific Liquid Water Content with Pressure'][16:37] #kg/kg
ecmwf_tciw_plevel_g = epg['Specific Ice Water Content with Pressure'][10:37] #kg/kg

#---ECMWF Southern Ocean Profile---#

ecmwf_tcc_alt_so = epso['Cloud Fraction'][11:37] # 0-1
ecmwf_tclw_alt_so = epso['Specific Liquid Water Content'][18:37] #kg/kg
ecmwf_tciw_alt_so = epso['Specific Ice Water Content'][11:37] #kg/kg
ecmwf_temp_alt_so = epso['Temperature Profile'][:] #K
ecmwf_plevel_alt_so = epso['Pressure Profile'][:] #hPa

ecmwf_tcc_temp_so = epso['Cloud Fraction with Temperature'][11:37] # 0-1
ecmwf_tclw_temp_so = epso['Specific Liquid Water Content with Temperature'][18:37] #kg/kg
ecmwf_tciw_temp_so = epso['Specific Ice Water Content with Temperature'][11:37] #kg/kg

ecmwf_tcc_plevel_so = epso['Cloud Fraction with Pressure'][12:37] # 0-1
ecmwf_tclw_plevel_so = epso['Specific Liquid Water Content with Pressure'][18:37] #kg/kg
ecmwf_tciw_plevel_so = epso['Specific Ice Water Content with Pressure'][11:37] #kg/kg

############################################################################### CCCM Data

#---CCCM Global Latitude Data---#

cccm_tcc_lat_g = c['tcc'][:] # 0-1
cccm_tclw_lat_g = c['tclw'][:] #kg/kg
cccm_tciw_lat_g = c['tciw'][:] #kg/kg

#---CCCM Southern Ocean Latitude Data---#

cccm_tcc_lat_so = cccm_tcc_lat_g[cccm_tcc_lat_g[:,0]>=-70]
cccm_tcc_lat_so = cccm_tcc_lat_so[cccm_tcc_lat_so[:,0]<=-50] # 0-1

cccm_tclw_lat_so = cccm_tclw_lat_g[cccm_tclw_lat_g[:,0]>=-70]
cccm_tclw_lat_so = cccm_tclw_lat_so[cccm_tclw_lat_so[:,0]<=-50] # kg/kg

cccm_tciw_lat_so = cccm_tciw_lat_g[cccm_tciw_lat_g[:,0]>=-70]
cccm_tciw_lat_so = cccm_tciw_lat_so[cccm_tciw_lat_so[:,0]<=-50] # kg/kg

#---CCCM Global Profile---#

cccm_tcc_alt_g = c['cf'][12:109] # 0-1
cccm_tclw_alt_g = c['lw'][70:133] #kg/kg
cccm_tciw_alt_g = c['iw'][36:133] #kg/kg
cccm_temp_alt_g = c['temp'][:] #K
cccm_plevel_alt_g = c['pressure'][:] #hPa

cccm_tcc_temp_g = c['cf_t'][:] # 0-1
cccm_tclw_temp_g = c['lw_t'][:] #kg/kg
cccm_tciw_temp_g = c['iw_t'][:] #kg/kg

#cccm_tcc_plevel_g = c['Cloud Fraction with Pressure'][12:110] # 0-1
#cccm_tclw_plevel_g = c['Specific Liquid Water Content with Pressure'][27:134] #kg/kg
#cccm_tciw_plevel_g = c['Specific Ice Water Content with Pressure'][36:134] #kg/kg

#---CCCM Southern Ocean Profile---#

cccm_tcc_alt_so = c['cf_so'][:] # 0-1
cccm_tclw_alt_so = c['lw_so'][53:] #kg/kg
cccm_tciw_alt_so = c['iw_so'][22:] #kg/kg
cccm_temp_alt_so = c['temp_so'][:] #K
cccm_plevel_alt_so = c['pressure_so'][:] #hPa

cccm_tcc_temp_so = c['cf_t_so'][:] # 0-1
cccm_tclw_temp_so = c['lw_t_so'][:] #kg/kg
cccm_tciw_temp_so = c['iw_t_so'][:] #kg/kg

#cccm_tcc_plevel_so = c['Cloud Fraction with Pressure'][:] # 0-1
#cccm_tclw_plevel_so = c['Specific Liquid Water Content with Pressure'][:] #kg/kg
#cccm_tciw_plevel_so = c['Specific Ice Water Content with Pressure'][:] #kg/kg

############################################################################### gfdl Data

#---gfdl Global Latitude Data---#

gfdl_tcc_lat_g = glg['Cloud Fraction'][:] # 0-1
gfdl_tclw_lat_g = glg['Liquid Water Content'][:] #kg/kg
gfdl_tciw_lat_g = glg['Ice Water Content'][:] #kg/kg

#---gfdl Southern Ocean Latitude Data---#

gfdl_tcc_lat_so = gfdl_tcc_lat_g[gfdl_tcc_lat_g[:,0]>=-70]
gfdl_tcc_lat_so = gfdl_tcc_lat_so[gfdl_tcc_lat_so[:,0]<=-50] # 0-1

gfdl_tclw_lat_so = gfdl_tclw_lat_g[gfdl_tclw_lat_g[:,0]>=-70]
gfdl_tclw_lat_so = gfdl_tclw_lat_so[gfdl_tclw_lat_so[:,0]<=-50] # kg/kg

gfdl_tciw_lat_so = gfdl_tciw_lat_g[gfdl_tciw_lat_g[:,0]>=-70]
gfdl_tciw_lat_so = gfdl_tciw_lat_so[gfdl_tciw_lat_so[:,0]<=-50] # kg/kg

#---gfdl Global Profile---#

gfdl_tcc_alt_g = gpg['Cloud Fraction'][:26] # 0-1
gfdl_tclw_alt_g = gpg['Specific Liquid Water Content'][:21] #kg/kg
gfdl_tciw_alt_g = gpg['Specific Ice Water Content'][:25] #kg/kg
gfdl_temp_alt_g = gpg['Temperature Profile'][:] #K
gfdl_plevel_alt_g = gpg['Pressure Profile'][:] #hPa

gfdl_tcc_temp_g = gpg['Cloud Fraction with Temperature'][:] # 0-1
gfdl_tclw_temp_g = gpg['Specific Liquid Water Content with Temperature'][:21] #kg/kg
gfdl_tciw_temp_g = gpg['Specific Ice Water Content with Temperature'][:] #kg/kg

gfdl_tcc_plevel_g = gpg['Cloud Fraction with Pressure'][:] # 0-1
gfdl_tclw_plevel_g = gpg['Specific Liquid Water Content with Pressure'][:21] #kg/kg
gfdl_tciw_plevel_g = gpg['Specific Ice Water Content with Pressure'][:] #kg/kg

#---gfdl Southern Ocean Profile---#

gfdl_tcc_alt_so = gpso['Cloud Fraction'][:23] # 0-1
gfdl_tclw_alt_so = gpso['Specific Liquid Water Content'][:19] #kg/kg
gfdl_tciw_alt_so = gpso['Specific Ice Water Content'][:23] #kg/kg
gfdl_temp_alt_so = gpso['Temperature Profile'][:] #K
gfdl_plevel_alt_so = gpso['Pressure Profile'][:] #hPa

gfdl_tcc_temp_so = gpso['Cloud Fraction with Temperature'][:23] # 0-1
gfdl_tclw_temp_so = gpso['Specific Liquid Water Content with Temperature'][:19] #kg/kg
gfdl_tciw_temp_so = gpso['Specific Ice Water Content with Temperature'][:23] #kg/kg

gfdl_tcc_plevel_so = gpso['Cloud Fraction with Pressure'][:23] # 0-1
gfdl_tclw_plevel_so = gpso['Specific Liquid Water Content with Pressure'][:19] #kg/kg
gfdl_tciw_plevel_so = gpso['Specific Ice Water Content with Pressure'][:23] #kg/kg



############################################################################### End importing Data

end = time.time()
print('Importing data took:', end - start, 's')

############################################################################### Global Latitude Plots

#---Plot Global Cloud Fraction with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tcc_lat_g[:,0],cccm_tcc_lat_g[:,1], '-r', label='CCCM')
ax.plot(ecmwf_tcc_lat_g[:,0],ecmwf_tcc_lat_g[:,1], '-b', label='ECMWF')
ax.plot(gfdl_tcc_lat_g[:,0],gfdl_tcc_lat_g[:,1], '-g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Fraction')

plt.title('2010 Global Cloud Fraction vs Latitude')

plt.grid(True)
plt.savefig("2010_tcc_lat_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Global Specific Liquid Water Content with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_lat_g[:,0],cccm_tclw_lat_g[:,1], '-r', label='CCCM')
ax.plot(ecmwf_tclw_lat_g[:,0],ecmwf_tclw_lat_g[:,1], '-b', label='ECMWF')
ax.plot(gfdl_tclw_lat_g[:,0],gfdl_tclw_lat_g[:,1], '-g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Content $kgm^{-2}$')

plt.title('2010 Global Liquid Water Content vs Latitude')

plt.grid(True)
plt.savefig("2010_tclw_lat_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Specific Ice Water Content with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tciw_lat_g[:,0],cccm_tciw_lat_g[:,1], '--r', label='CCCM')
ax.plot(ecmwf_tciw_lat_g[:,0],ecmwf_tciw_lat_g[:,1], '--b', label='ECMWF')
ax.plot(gfdl_tciw_lat_g[:,0],gfdl_tciw_lat_g[:,1], '--g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Content $kgm^{-2}$')

plt.title('2010 Global Ice Water Content vs Latitude')

plt.grid(True)
plt.savefig("2010_tciw_lat_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Comparison Specific Liquid and Ice Water Content with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_lat_g[:,0],cccm_tclw_lat_g[:,1], '-r', label='Liquid - CCCM')
ax.plot(ecmwf_tclw_lat_g[:,0],ecmwf_tclw_lat_g[:,1], '-b', label='Liquid - ECMWF')
ax.plot(gfdl_tclw_lat_g[:,0],gfdl_tclw_lat_g[:,1], '-g', label='Liquid - GFDL')
ax.plot(cccm_tciw_lat_g[:,0],cccm_tciw_lat_g[:,1], '--r', label='Ice - CCCM')
ax.plot(ecmwf_tciw_lat_g[:,0],ecmwf_tciw_lat_g[:,1], '--b', label='Ice - ECMWF')
ax.plot(gfdl_tciw_lat_g[:,0],gfdl_tciw_lat_g[:,1], '--g', label='Ice - GFDL')

ax.legend(loc='upper center', bbox_to_anchor=(1.2, 1.0));

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Content $kgm^{-2}$')

plt.title('2010 Global Liquid and Ice Water Content vs Latitude')

plt.grid(True)
plt.savefig("2010_tclw_tciw_lat_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CCCM Cloud Fraction and Phase Content with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(cccm_tcc_lat_g[:,0],cccm_tcc_lat_g[:,1], '-r', label='Cloud Fraction')
ax2.plot(cccm_tclw_lat_g[:,0],cccm_tclw_lat_g[:,1], '-b', label='Liquid Content')
ax2.plot(cccm_tciw_lat_g[:,0],cccm_tciw_lat_g[:,1], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 CCCM Global Cloud Fraction and Phase Content vs Latitude')

plt.grid(True)
plt.savefig("2010_CCCM_lat_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot ECMWF Cloud Fraction and Phase Content with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(ecmwf_tcc_lat_g[:,0],ecmwf_tcc_lat_g[:,1], '-r', label='Cloud Fraction')
ax2.plot(ecmwf_tclw_lat_g[:,0],ecmwf_tclw_lat_g[:,1]*10000, '-b', label='Liquid Content')
ax2.plot(ecmwf_tciw_lat_g[:,0],ecmwf_tciw_lat_g[:,1]*10000, '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 ECMWF Global Cloud Fraction and Phase Content vs Latitude')

plt.grid(True)
plt.savefig("2010_ECMWF_lat_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""
#---Plot GFDL Cloud Fraction and Phase Content with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(gfdl_tcc_lat_g[:,0],gfdl_tcc_lat_g[:,1], '-r', label='Cloud Fraction')
ax2.plot(gfdl_tclw_lat_g[:,0],gfdl_tclw_lat_g[:,1]*10000, '-b', label='Liquid Content')
ax2.plot(gfdl_tciw_lat_g[:,0],gfdl_tciw_lat_g[:,1]*10000, '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 GFDL Global Cloud Fraction and Phase Content vs Latitude')

plt.grid(True)
plt.savefig("2010_GFDL_lat_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""
############################################################################### Southern Ocean Latitude Plots

#---Plot Southern Ocean Cloud Fraction with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tcc_lat_so[:,0],cccm_tcc_lat_so[:,1], '-r', label='CCCM')
ax.plot(ecmwf_tcc_lat_so[:,0],ecmwf_tcc_lat_so[:,1], '-b', label='ECMWF')
ax.plot(gfdl_tcc_lat_so[:,0],gfdl_tcc_lat_so[:,1], '-g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Fraction')

plt.title('2010 Southern Ocean Cloud Fraction vs Latitude')

plt.grid(True)
plt.savefig("2010_tcc_lat_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Southern Ocean Specific Liquid Water Content with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_lat_so[:,0],cccm_tclw_lat_so[:,1]*10000, '-r', label='CCCM')
ax.plot(ecmwf_tclw_lat_so[:,0],ecmwf_tclw_lat_so[:,1]*10000, '-b', label='ECMWF')
ax.plot(gfdl_tclw_lat_so[:,0],gfdl_tclw_lat_so[:,1]*10000, '-g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Latitude')
ax.set_ylabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Southern Ocean Specific Liquid Water Content vs Latitude')

plt.grid(True)
plt.savefig("2010_tclw_lat_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Southern Ocean Specific Ice Water Content with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tciw_lat_so[:,0],cccm_tciw_lat_so[:,1]*10000, '--r', label='CCCM')
ax.plot(ecmwf_tciw_lat_so[:,0],ecmwf_tciw_lat_so[:,1]*10000, '--b', label='ECMWF')
ax.plot(gfdl_tciw_lat_so[:,0],gfdl_tciw_lat_so[:,1]*10000, '--g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Latitude')
ax.set_ylabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Southern Ocean Specific Ice Water Content vs Latitude')

plt.grid(True)
plt.savefig("2010_tciw_lat_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Comparison Specific Liquid and Ice Water Content with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_lat_so[:,0],cccm_tclw_lat_so[:,1]*10000, '-r', label='Liquid - CCCM')
ax.plot(ecmwf_tclw_lat_so[:,0],ecmwf_tclw_lat_so[:,1]*10000, '-b', label='Liquid - ECMWF')
ax.plot(gfdl_tclw_lat_so[:,0],gfdl_tclw_lat_so[:,1]*10000, '-g', label='Liquid - GFDL')
ax.plot(cccm_tciw_lat_so[:,0],cccm_tciw_lat_so[:,1]*10000, '--r', label='Ice - CCCM')
ax.plot(ecmwf_tciw_lat_so[:,0],ecmwf_tciw_lat_so[:,1]*10000, '--b', label='Ice - ECMWF')
ax.plot(gfdl_tciw_lat_so[:,0],gfdl_tciw_lat_so[:,1]*10000, '--g', label='Ice - GFDL')

ax.legend(loc='upper center', bbox_to_anchor=(1.2, 1.0));

ax.set_xlabel('Latitude')
ax.set_ylabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Southern Ocean Specific Liquid and Ice Water Content vs Latitude')

plt.grid(True)

plt.savefig("2010_tclw_tciw_lat_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""
#---Plot CCCM Cloud Fraction and Phase Content with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(cccm_tcc_lat_so[:,0],cccm_tcc_lat_so[:,1], '-r', label='Cloud Fraction')
ax2.plot(cccm_tclw_lat_so[:,0],cccm_tclw_lat_so[:,1], '-b', label='Liquid Content')
ax2.plot(cccm_tciw_lat_so[:,0],cccm_tciw_lat_so[:,1], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 CCCM Southern Ocean Cloud Fraction and Phase Content vs Latitude')


plt.savefig("2010_CCCM_lat_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot ECMWF Cloud Fraction and Phase Content with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(ecmwf_tcc_lat_so[:,0],ecmwf_tcc_lat_so[:,1], '-r', label='Cloud Fraction')
ax2.plot(ecmwf_tclw_lat_so[:,0],ecmwf_tclw_lat_so[:,1]*10000, '-b', label='Liquid Content')
ax2.plot(ecmwf_tciw_lat_so[:,0],ecmwf_tciw_lat_so[:,1]*10000, '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 ECMWF Southern Ocean Cloud Fraction and Phase Content vs Latitude')


plt.savefig("2010_ECMWF_lat_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""
#---Plot GFDL Cloud Fraction and Phase Content with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(gfdl_tcc_lat_so[:,0],gfdl_tcc_lat_so[:,1], '-r', label='Cloud Fraction')
ax2.plot(gfdl_tclw_lat_so[:,0],gfdl_tclw_lat_so[:,1]*10000, '-b', label='Liquid Content')
ax2.plot(gfdl_tciw_lat_so[:,0],gfdl_tciw_lat_so[:,1]*10000, '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 GFDL Southern Ocean Cloud Fraction and Phase Content vs Latitude')

plt.savefig("2010_GFDL_lat_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""



############################################################################### Global Altitude Plots

#---Plot Global Cloud Fraction Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tcc_alt_g[:,1],cccm_tcc_alt_g[:,0], '-r', label='CCCM')
ax.plot(ecmwf_tcc_alt_g[:,1],ecmwf_tcc_alt_g[:,0], '-b', label='ECMWF')
ax.plot(gfdl_tcc_alt_g[:,1],gfdl_tcc_alt_g[:,0], '-g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Cloud Fraction')

plt.title('2010 Global Cloud Fraction vs Altitude')

plt.grid(True)

plt.savefig("2010_tcc_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Specific Liquid Water Content Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_alt_g[:,1]*10000,cccm_tclw_alt_g[:,0], '-r', label='CCCM')
ax.plot(ecmwf_tclw_alt_g[:,1]*10000,ecmwf_tclw_alt_g[:,0], '-b', label='ECMWF')
ax.plot(gfdl_tclw_alt_g[:,1]*10000,gfdl_tclw_alt_g[:,0], '-g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Global Specific Liquid Water Content vs Altitude')

plt.grid(True)
plt.savefig("2010_tclw_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Specific Ice Water Content Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tciw_alt_g[:,1]*10000,cccm_tciw_alt_g[:,0], '--r', label='CCCM')
ax.plot(ecmwf_tciw_alt_g[:,1]*10000,ecmwf_tciw_alt_g[:,0], '--b', label='ECMWF')
ax.plot(gfdl_tciw_alt_g[:,1]*10000,gfdl_tciw_alt_g[:,0], '--g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Global Specific Ice Water Content vs Altitude')

plt.grid(True)
plt.savefig("2010_tciw_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Comparison Specific Liquid and Ice Water Content Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_alt_g[:,1]*10000,cccm_tclw_alt_g[:,0], '-r', label='Liquid - CCCM')
ax.plot(ecmwf_tclw_alt_g[:,1]*10000,ecmwf_tclw_alt_g[:,0], '-b', label='Liquid - ECMWF')
ax.plot(gfdl_tclw_alt_g[:,1]*10000,gfdl_tclw_alt_g[:,0], '-g', label='Liquid - GFDL')
ax.plot(cccm_tciw_alt_g[:,1]*10000,cccm_tciw_alt_g[:,0], '--r', label='Ice - CCCM')
ax.plot(ecmwf_tciw_alt_g[:,1]*10000,ecmwf_tciw_alt_g[:,0], '--b', label='Ice - ECMWF')
ax.plot(gfdl_tciw_alt_g[:,1]*10000,gfdl_tciw_alt_g[:,0], '--g', label='Ice - GFDL')

ax.legend(loc='upper center', bbox_to_anchor=(1.2, 1.0));

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Global Specific Liquid and Ice Water Content vs Altitude')

plt.grid(True)
plt.savefig("2010_tclw_tciw_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Temperature Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_temp_alt_g[:,1],cccm_temp_alt_g[:,0], '-r', label='CCCM')
ax.plot(ecmwf_temp_alt_g[:,1],ecmwf_temp_alt_g[:,0], '-b', label='ECMWF')
ax.plot(gfdl_temp_alt_g[:,1],gfdl_temp_alt_g[:,0], '-g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Temperature (K)')

plt.title('2010 Global Temperature vs Altitude')

plt.grid(True)
plt.savefig("2010_T_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Pressure Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_plevel_alt_g[:,1],cccm_plevel_alt_g[:,0], '-r', label='CCCM')
ax.plot(ecmwf_plevel_alt_g[:,1],ecmwf_plevel_alt_g[:,0], '-b', label='ECMWF')
ax.plot(gfdl_plevel_alt_g[:,1],gfdl_plevel_alt_g[:,0], '-g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Pressure (hPa)')

plt.title('2010 Global Temperature vs Altitude')
#plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("2010_P_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Specific Liquid Water Content with Temperature Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_temp_g[:,1]*10000,cccm_tclw_temp_g[:,0], '-r', label='CCCM')
ax.plot(ecmwf_tclw_temp_g[:,1]*10000,ecmwf_tclw_temp_g[:,0], '-b', label='ECMWF')
ax.plot(gfdl_tclw_temp_g[:,1]*10000,gfdl_tclw_temp_g[:,0], '-g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Temperature (K)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Global Specific Liquid Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("2010_tclw_T_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Specific Ice Water Content with Temperature Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tciw_temp_g[:,1]*10000,cccm_tciw_temp_g[:,0], '--r', label='CCCM')
ax.plot(ecmwf_tciw_temp_g[:,1]*10000,ecmwf_tciw_temp_g[:,0], '--b', label='ECMWF')
ax.plot(gfdl_tciw_temp_g[:,1]*10000,gfdl_tciw_temp_g[:,0], '--g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Temperature (K)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Global Specific Ice Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("2010_tciw_T_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Comparison Specific Liquid and Ice Water Content with Temperature Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_temp_g[:,1]*10000,cccm_tclw_temp_g[:,0], '-r', label='Liquid - CCCM')
ax.plot(ecmwf_tclw_temp_g[:,1]*10000,ecmwf_tclw_temp_g[:,0], '-b', label='Liquid - ECMWF')
ax.plot(gfdl_tclw_temp_g[:,1]*10000,gfdl_tclw_temp_g[:,0], '-g', label='Liquid - GFDL')
ax.plot(cccm_tciw_temp_g[:,1]*10000,cccm_tciw_temp_g[:,0], '--r', label='Ice - CCCM')
ax.plot(ecmwf_tciw_temp_g[:,1]*10000,ecmwf_tciw_temp_g[:,0], '--b', label='Ice - ECMWF')
ax.plot(gfdl_tciw_temp_g[:,1]*10000,gfdl_tciw_temp_g[:,0], '--g', label='Ice - GFDL')

ax.legend(loc='upper center', bbox_to_anchor=(1.2, 1.0));

ax.set_ylabel('Temperature (K)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Global Specific Liquid and Ice Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("2010_tclw_tciw_T_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""
#---Plot CCCM Cloud Fraction and Phase Profile---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()
ax1.plot(cccm_tcc_alt_g[:,1],cccm_tcc_alt_g[:,0], '-r', label='Cloud Fraction')
ax2.plot(cccm_tclw_alt_g[:,1]*10000,cccm_tclw_alt_g[:,0], '-b', label='Liquid Content')
ax2.plot(cccm_tciw_alt_g[:,1]*10000,cccm_tciw_alt_g[:,0], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_ylabel('Altitude')
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 CCCM Global Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("2010_CCCM_alt_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""
#---Plot ECMWF Cloud Fraction and Phase Profile---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()
ax1.plot(ecmwf_tcc_alt_g[:,1],ecmwf_tcc_alt_g[:,0], '-r', label='Cloud Fraction')
ax2.plot(ecmwf_tclw_alt_g[:,1]*10000,ecmwf_tclw_alt_g[:,0], '-b', label='Liquid Content')
ax2.plot(ecmwf_tciw_alt_g[:,1]*10000,ecmwf_tciw_alt_g[:,0], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_ylabel('Altitude')
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 ECMWF Global Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("2010_ECMWF_alt_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot GFDL Cloud Fraction and Phase Profile---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()
ax1.plot(gfdl_tcc_alt_g[:,1],gfdl_tcc_alt_g[:,0], '-r', label='Cloud Fraction')
ax2.plot(gfdl_tclw_alt_g[:,1]*10000,gfdl_tclw_alt_g[:,0], '-b', label='Liquid Content')
ax2.plot(gfdl_tciw_alt_g[:,1]*10000,gfdl_tciw_alt_g[:,0], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_ylabel('Altitude')
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 GFDL Global Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("2010_GFDL_alt_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

############################################################################### Southern Ocean Altitude Plots

#---Plot Southern Ocean Cloud Fraction Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tcc_alt_so[:,1],cccm_tcc_alt_so[:,0], '-r', label='CCCM')
ax.plot(ecmwf_tcc_alt_so[:,1],ecmwf_tcc_alt_so[:,0], '-b', label='ECMWF')
ax.plot(gfdl_tcc_alt_so[:,1],gfdl_tcc_alt_so[:,0], '-g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Cloud Fraction')

plt.title('2010 Southern Ocean Cloud Fraction vs Altitude')

plt.grid(True)
plt.savefig("2010_tcc_alt_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Specific Liquid Water Content Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_alt_so[:,1]*10000,cccm_tclw_alt_so[:,0], '-r', label='CCCM')
ax.plot(ecmwf_tclw_alt_so[:,1]*10000,ecmwf_tclw_alt_so[:,0], '-b', label='ECMWF')
ax.plot(gfdl_tclw_alt_so[:,1]*10000,gfdl_tclw_alt_so[:,0], '-g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Southern Ocean Specific Liquid Water Content vs Altitude')

plt.grid(True)
plt.savefig("2010_tclw_alt_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Specific Ice Water Content Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tciw_alt_so[:,1]*10000,cccm_tciw_alt_so[:,0], '--r', label='CCCM')
ax.plot(ecmwf_tciw_alt_so[:,1]*10000,ecmwf_tciw_alt_so[:,0], '--b', label='ECMWF')
ax.plot(gfdl_tciw_alt_so[:,1]*10000,gfdl_tciw_alt_so[:,0], '--g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Southern Ocean Specific Ice Water Content vs Altitude')

plt.grid(True)
plt.savefig("2010_tciw_alt_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Comparison Specific Liquid and Ice Water Content Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_alt_so[:,1]*10000,cccm_tclw_alt_so[:,0], '-r', label='Liquid - CCCM')
ax.plot(ecmwf_tclw_alt_so[:,1]*10000,ecmwf_tclw_alt_so[:,0], '-b', label='Liquid - ECMWF')
ax.plot(gfdl_tclw_alt_so[:,1]*10000,gfdl_tclw_alt_so[:,0], '-g', label='Liquid - GFDL')
ax.plot(cccm_tciw_alt_so[:,1]*10000,cccm_tciw_alt_so[:,0], '--r', label='Ice - CCCM')
ax.plot(ecmwf_tciw_alt_so[:,1]*10000,ecmwf_tciw_alt_so[:,0], '--b', label='Ice - ECMWF')
ax.plot(gfdl_tciw_alt_so[:,1]*10000,gfdl_tciw_alt_so[:,0], '--g', label='Ice - GFDL')

ax.legend(loc='upper center', bbox_to_anchor=(1.2, 1.0));

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Southern Ocean Specific Liquid and Ice Water Content vs Altitude')

plt.grid(True)
plt.savefig("2010_tclw_tciw_alt_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Temperature Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_temp_alt_so[:,1],cccm_temp_alt_so[:,0], '-r', label='CCCM')
ax.plot(ecmwf_temp_alt_so[:,1],ecmwf_temp_alt_so[:,0], '-b', label='ECMWF')
ax.plot(gfdl_temp_alt_so[:,1],gfdl_temp_alt_so[:,0], '-g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Temperature (K)')

plt.title('2010 Southern Ocean Temperature vs Altitude')

plt.grid(True)
plt.savefig("2010_T_alt_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Pressure Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_plevel_alt_so[:,1],cccm_plevel_alt_so[:,0], '-r', label='CCCM')
ax.plot(ecmwf_plevel_alt_so[:,1],ecmwf_plevel_alt_so[:,0], '-b', label='ECMWF')
ax.plot(gfdl_plevel_alt_so[:,1],gfdl_plevel_alt_so[:,0], '-g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Pressure (hPa)')

plt.title('2010 Southern Ocean Temperature vs Altitude')
#plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("2010_P_alt_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Specific Liquid Water Content with Temperature Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_temp_so[:,1]*10000,cccm_tclw_temp_so[:,0], '-r', label='CCCM')
ax.plot(ecmwf_tclw_temp_so[:,1]*10000,ecmwf_tclw_temp_so[:,0], '-b', label='ECMWF')
ax.plot(gfdl_tclw_temp_so[:,1]*10000,gfdl_tclw_temp_so[:,0], '-g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Temperature (K)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Southern Ocean Specific Liquid Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("2010_tclw_T_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Specific Ice Water Content with Temperature Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tciw_temp_so[:,1]*10000,cccm_tciw_temp_so[:,0], '--r', label='CCCM')
ax.plot(ecmwf_tciw_temp_so[:,1]*10000,ecmwf_tciw_temp_so[:,0], '--b', label='ECMWF')
ax.plot(gfdl_tciw_temp_so[:,1]*10000,gfdl_tciw_temp_so[:,0], '--g', label='GFDL')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Temperature (K)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Southern Ocean Specific Ice Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("2010_tciw_T_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Comparison Specific Liquid and Ice Water Content with Temperature Profile---#

"""
plt.figure()
fig, ax = plt.subplots()


ax.plot(cccm_tclw_temp_so[:,1]*10000,cccm_tclw_temp_so[:,0], '-r', label='Liquid - CCCM')
ax.plot(ecmwf_tclw_temp_so[:,1]*10000,ecmwf_tclw_temp_so[:,0], '-b', label='Liquid - ECMWF')
ax.plot(gfdl_tclw_temp_so[:,1]*10000,gfdl_tclw_temp_so[:,0], '-g', label='Liquid - GFDL')
ax.plot(cccm_tciw_temp_so[:,1]*10000,cccm_tciw_temp_so[:,0], '--r', label='Ice - CCCM')
ax.plot(ecmwf_tciw_temp_so[:,1]*10000,ecmwf_tciw_temp_so[:,0], '--b', label='Ice - ECMWF')
ax.plot(gfdl_tciw_temp_so[:,1]*10000,gfdl_tciw_temp_so[:,0], '--g', label='Ice - GFDL')

ax.legend(loc='upper center', bbox_to_anchor=(1.2, 1.0));

ax.set_ylabel('Temperature (K)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Southern Ocean Specific Liquid and Ice Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("2010_tclw_tciw_T_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CCCM Cloud Fraction and Phase Profile---#
"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()
ax1.plot(cccm_tcc_alt_so[:,1],cccm_tcc_alt_so[:,0], '-r', label='Cloud Fraction')
ax2.plot(cccm_tclw_alt_so[:,1]*10000,cccm_tclw_alt_so[:,0], '-b', label='Liquid Content')
ax2.plot(cccm_tciw_alt_so[:,1]*10000,cccm_tciw_alt_so[:,0], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_ylabel('Altitude')
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 CCCM Southern Ocean Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("2010_CCCM_alt_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""
#---Plot ECMWF Cloud Fraction and Phase Profile---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()
ax1.plot(ecmwf_tcc_alt_so[:,1],ecmwf_tcc_alt_so[:,0], '-r', label='Cloud Fraction')
ax2.plot(ecmwf_tclw_alt_so[:,1]*10000,ecmwf_tclw_alt_so[:,0], '-b', label='Liquid Content')
ax2.plot(ecmwf_tciw_alt_so[:,1]*10000,ecmwf_tciw_alt_so[:,0], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_ylabel('Altitude')
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 ECMWF Southern Ocean Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("2010_ECMWF_alt_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot GFDL Cloud Fraction and Phase Profile---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()
ax1.plot(gfdl_tcc_alt_so[:,1],gfdl_tcc_alt_so[:,0], '-r', label='Cloud Fraction')
ax2.plot(gfdl_tclw_alt_so[:,1]*10000,gfdl_tclw_alt_so[:,0], '-b', label='Liquid Content')
ax2.plot(gfdl_tciw_alt_so[:,1]*10000,gfdl_tciw_alt_so[:,0], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_ylabel('Altitude')
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('2010 GFDL Southern Ocean Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("2010_GFDL_alt_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""


#---Plot Global and Southern Ocean Comparison Specific Liquid and Ice Water Content vs Altitude---#

"""
plt.figure()
fig, ax = plt.subplots()

#ax.plot(cccm_tclw_alt_g[:,1],cccm_tclw_alt_g[:,0], '-r', label='Liquid - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tclw_alt_g[:,1],ecmwf_tclw_alt_g[:,0], '-b', label='Liquid - ECMWF.ERA5 Reanalysis Model')
#ax.plot(gfdl_tclw_alt_g[:,1],gfdl_tclw_alt_g[:,0], '-g', label='Liquid - GFDL')
#ax.plot(cccm_tciw_alt_g[:,1],cccm_tciw_alt_g[:,0], '--r', label='Ice - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tciw_alt_g[:,1],ecmwf_tciw_alt_g[:,0], '--b', label='Ice - ECMWF.ERA5 Reanalysis Model')
#ax.plot(gfdl_tciw_alt_g[:,1],gfdl_tciw_alt_g[:,0], '--g', label='Ice - GFDL')
#ax.plot(cccm_tclw_alt_so[:,1],cccm_tclw_alt_so[:,0], '-r', label='Liquid - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tclw_alt_so[:,1],ecmwf_tclw_alt_so[:,0], '-b', label='Liquid - ECMWF.ERA5 Reanalysis Model')
#ax.plot(gfdl_tclw_alt_so[:,1],gfdl_tclw_alt_so[:,0], '-g', label='Liquid - GFDL')
#ax.plot(cccm_tciw_alt_so[:,1],cccm_tciw_alt_so[:,0], '--r', label='Ice - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tciw_alt_so[:,1],ecmwf_tciw_alt_so[:,0], '--b', label='Ice - ECMWF.ERA5 Reanalysis Model')
#ax.plot(gfdl_tciw_alt_so[:,1],gfdl_tciw_alt_so[:,0], '--g', label='Ice - GFDL')

ax.legend(loc='upper center', bbox_to_anchor=(1.2, 1.0));

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Specific Liquid and Ice Water Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Southern Ocean Specific Liquid and Ice Water Content vs Altitude')

plt.grid(True)
plt.show()
"""
#---Plot Global and Southern Ocean Comparison Specific Liquid and Ice Water Content with Temperature Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

#ax.plot(cccm_tclw_temp_g[:,1],cccm_tclw_temp_g[:,0], '-r', label='Liquid - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tclw_temp_g[:,1],ecmwf_tclw_temp_g[:,0], '-b', label='Liquid - ECMWF.ERA5 Reanalysis Model')
#ax.plot(gfdl_tclw_temp_g[:,1],gfdl_tclw_temp_g[:,0], '-g', label='Liquid - GFDL')
#ax.plot(cccm_tciw_temp_g[:,1],cccm_tciw_temp_g[:,0], '--r', label='Ice - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tciw_temp_g[:,1],ecmwf_tciw_temp_g[:,0], '--b', label='Ice - ECMWF.ERA5 Reanalysis Model')
#ax.plot(gfdl_tciw_temp_g[:,1],gfdl_tciw_temp_g[:,0], '--g', label='Ice - GFDL')
#ax.plot(cccm_tclw_temp_so[:,1],cccm_tclw_temp_so[:,0], '-r', label='Liquid - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tclw_temp_so[:,1],ecmwf_tclw_temp_so[:,0], '-b', label='Liquid - ECMWF.ERA5 Reanalysis Model')
#ax.plot(gfdl_tclw_temp_so[:,1],gfdl_tclw_temp_so[:,0], '-g', label='Liquid - GFDL')
#ax.plot(cccm_tciw_temp_so[:,1],cccm_tciw_temp_so[:,0], '--r', label='Ice - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tciw_temp_so[:,1],ecmwf_tciw_temp_so[:,0], '--b', label='Ice - ECMWF.ERA5 Reanalysis Model')
#ax.plot(gfdl_tciw_temp_so[:,1],gfdl_tciw_temp_so[:,0], '--g', label='Ice - GFDL')

ax.legend(loc='upper center', bbox_to_anchor=(1.2, 1.0));

ax.set_ylabel('Temperature (K)')
ax.set_xlabel('Specific Liquid and Ice Water Content $(kg/kg) x 10^{-4}$')

plt.title('2010 Southern Ocean Specific Liquid and Ice Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.show()
"""