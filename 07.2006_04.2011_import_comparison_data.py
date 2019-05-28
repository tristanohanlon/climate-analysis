# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

This import the reduced datasets from 07.2006 to 04.2011 CCCM, ECMWF and GFDL. 
The code can select both global or southern ocean data.

"""
import time
import matplotlib.pyplot as plt
import h5py
import os
import numpy as np

start = time.time()

"""
#---Importing Data from Reduced Datasets---#

# Uni Laptop
#ECMWF Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets')
a = h5py.File('07.2006_04.2011_ECMWF.h5', 'r')

#CCCM Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
c = h5py.File('07.2006_04.2011_CCCM.h5', 'r')

#GFDL Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GFDL/reduced_datasets')
b = h5py.File('07.2006_04.2011_gfdl.h5', 'r')

#CAM5 Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CAM5/reduced_datasets')
d = h5py.File('07.2006_04.2011_CAM5.h5', 'r')

"""

# Home PC
#ECMWF Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets')
a = h5py.File('07.2006_04.2011_ECMWF.h5', 'r')

#CCCM Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
c = h5py.File('07.2006_04.2011_CCCM.h5', 'r')

#GFDL Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL/reduced_datasets')
b = h5py.File('07.2006_04.2011_gfdl.h5', 'r')

#CAM5 Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CAM5/reduced_datasets')
d = h5py.File('07.2006_04.2011_CAM5.h5', 'r')

"""
# Laptop
#ECMWF Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets')
a = h5py.File('07.2006_04.2011_ECMWF.h5', 'r')

#CCCM Data
#os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
c = h5py.File('07.2006_04.2011_CCCM.h5', 'r')

#GFDL Data
#os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/GFDL/reduced_datasets')
b = h5py.File('07.2006_04.2011_gfdl.h5', 'r')

#CAM5 Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CAM5/reduced_datasets')
d = h5py.File('07.2006_04.2011_CAM5.h5', 'r')

"""

############################################################################### ECMWF Data

#---ECMWF Global Latitude Data---#

ecmwf_tcc_lat_g = a['tcc'][:] # 0-1
ecmwf_tclw_lat_g = a['tclw'][:] #kgm^-2
ecmwf_tciw_lat_g = a['tciw'][:] #kgm^-2

#---ECMWF Southern Ocean Latitude Data---#

ecmwf_tcc_lat_so = ecmwf_tcc_lat_g[ecmwf_tcc_lat_g[:,0]>=-70]
ecmwf_tcc_lat_so = ecmwf_tcc_lat_so[ecmwf_tcc_lat_so[:,0]<=-50] # 0-1

ecmwf_tclw_lat_so = ecmwf_tclw_lat_g[ecmwf_tclw_lat_g[:,0]>=-70]
ecmwf_tclw_lat_so = ecmwf_tclw_lat_so[ecmwf_tclw_lat_so[:,0]<=-50] #kgm^-2

ecmwf_tciw_lat_so = ecmwf_tciw_lat_g[ecmwf_tciw_lat_g[:,0]>=-70]
ecmwf_tciw_lat_so = ecmwf_tciw_lat_so[ecmwf_tciw_lat_so[:,0]<=-50] #kgm^-2

#---ECMWF Phase Fractions---#

ecmwf_tclw_frac_lat_g = (ecmwf_tclw_lat_g[:,1] / (ecmwf_tclw_lat_g[:,1] + ecmwf_tciw_lat_g[:,1])) * ecmwf_tcc_lat_g[:,1]
ecmwf_tciw_frac_lat_g = (ecmwf_tciw_lat_g[:,1] / (ecmwf_tclw_lat_g[:,1] + ecmwf_tciw_lat_g[:,1])) * ecmwf_tcc_lat_g[:,1]

ecmwf_tclw_frac_lat_g = np.vstack((ecmwf_tclw_lat_g[:,0], ecmwf_tclw_frac_lat_g)).T
ecmwf_tciw_frac_lat_g = np.vstack((ecmwf_tciw_lat_g[:,0], ecmwf_tciw_frac_lat_g)).T

#---ECMWF Global Profile---#

ecmwf_tcc_alt_g = a['cf'][8:] # 0-1
ecmwf_tclw_alt_g = a['lw'][14:] #kg/kg
ecmwf_tciw_alt_g = a['iw'][9:] #kg/kg
ecmwf_temp_alt_g = a['temp'][:] #K
ecmwf_plevel_alt_g = a['pressure'][:] #hPa

ecmwf_tcc_temp_g = a['cf_t'][13:] # 0-1
ecmwf_tclw_temp_g = a['lw_t'][17:] #kg/kg
ecmwf_tciw_temp_g = a['iw_t'][10:] #kg/kg

#ecmwf_tcc_plevel_g = a['Cloud Fraction with Pressure'][8:37] # 0-1
#ecmwf_tclw_plevel_g = a['Specific Liquid Water Content with Pressure'][16:37] #kg/kg
#ecmwf_tciw_plevel_g = a['Specific Ice Water Content with Pressure'][10:37] #kg/kg

#---ECMWF Southern Ocean Profile---#

ecmwf_tcc_alt_so = a['cf_so'][11:] # 0-1
ecmwf_tclw_alt_so = a['lw_so'][18:] #kg/kg
ecmwf_tciw_alt_so = a['iw_so'][11:] #kg/kg
ecmwf_temp_alt_so = a['temp_so'][:] #K
ecmwf_plevel_alt_so = a['pressure_so'][:] #hPa

ecmwf_tcc_temp_so = a['cf_t_so'][11:] # 0-1
ecmwf_tclw_temp_so = a['lw_t_so'][18:] #kg/kg
ecmwf_tciw_temp_so = a['iw_t_so'][11:] #kg/kg

#ecmwf_tcc_plevel_so = a['Cloud Fraction with Pressure'][12:37] # 0-1
#ecmwf_tclw_plevel_so = a['Specific Liquid Water Content with Pressure'][18:37] #kg/kg
#ecmwf_tciw_plevel_so = a['Specific Ice Water Content with Pressure'][11:37] #kg/kg

############################################################################### CCCM Data

#---CCCM Global Latitude Data---#

cccm_tcc_lat_g = c['tcc'][:] # 0-1
cccm_tclw_lat_g = c['tclw'][:] #kgm^-2
cccm_tciw_lat_g = c['tciw'][:] #kgm^-2

#---CCCM Southern Ocean Latitude Data---#

cccm_tcc_lat_so = cccm_tcc_lat_g[cccm_tcc_lat_g[:,0]>=-70]
cccm_tcc_lat_so = cccm_tcc_lat_so[cccm_tcc_lat_so[:,0]<=-50] # 0-1

cccm_tclw_lat_so = cccm_tclw_lat_g[cccm_tclw_lat_g[:,0]>=-70]
cccm_tclw_lat_so = cccm_tclw_lat_so[cccm_tclw_lat_so[:,0]<=-50] #kgm^-2

cccm_tciw_lat_so = cccm_tciw_lat_g[cccm_tciw_lat_g[:,0]>=-70]
cccm_tciw_lat_so = cccm_tciw_lat_so[cccm_tciw_lat_so[:,0]<=-50] #kgm^-2

#---CCCM Phase Fractions---#

cccm_tclw_frac_lat_g = (cccm_tclw_lat_g[:,1] / (cccm_tclw_lat_g[:,1] + cccm_tciw_lat_g[:,1])) * cccm_tcc_lat_g[:,1]
cccm_tciw_frac_lat_g = (cccm_tciw_lat_g[:,1] / (cccm_tclw_lat_g[:,1] + cccm_tciw_lat_g[:,1])) * cccm_tcc_lat_g[:,1]

cccm_tclw_frac_lat_g = np.vstack((cccm_tclw_lat_g[:,0], cccm_tclw_frac_lat_g)).T
cccm_tciw_frac_lat_g = np.vstack((cccm_tciw_lat_g[:,0], cccm_tciw_frac_lat_g)).T

#---CCCM Global Profile---#

cccm_tcc_alt_g = c['cf'][4:101] # 0-1
cccm_tclw_alt_g = c['lw'][4:63] #kg/kg
cccm_tciw_alt_g = c['iw'][4:93] #kg/kg
cccm_temp_alt_g = c['temp'][:] #K
cccm_plevel_alt_g = c['pressure'][:] #hPa

cccm_tcc_temp_g = c['cf_t'][:80] # 0-1
cccm_tclw_temp_g = c['lw_t'][:100] #kg/kg
cccm_tciw_temp_g = c['iw_t'][:100] #kg/kg

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

gfdl_tcc_lat_g = b['tcc'][:] # 0-1
gfdl_tclw_lat_g = b['tclw'][:] #kgm^-2
gfdl_tciw_lat_g = b['tciw'][:] #kgm^-2

#---gfdl Southern Ocean Latitude Data---#

gfdl_tcc_lat_so = gfdl_tcc_lat_g[gfdl_tcc_lat_g[:,0]>=-70]
gfdl_tcc_lat_so = gfdl_tcc_lat_so[gfdl_tcc_lat_so[:,0]<=-50] # 0-1

gfdl_tclw_lat_so = gfdl_tclw_lat_g[gfdl_tclw_lat_g[:,0]>=-70]
gfdl_tclw_lat_so = gfdl_tclw_lat_so[gfdl_tclw_lat_so[:,0]<=-50] #kgm^-2

gfdl_tciw_lat_so = gfdl_tciw_lat_g[gfdl_tciw_lat_g[:,0]>=-70]
gfdl_tciw_lat_so = gfdl_tciw_lat_so[gfdl_tciw_lat_so[:,0]<=-50] #kgm^-2

#---CCCM Phase Fractions---#

gfdl_tclw_frac_lat_g = (gfdl_tclw_lat_g[:,1] / (gfdl_tclw_lat_g[:,1] + gfdl_tciw_lat_g[:,1])) * gfdl_tcc_lat_g[:,1]
gfdl_tciw_frac_lat_g = (gfdl_tciw_lat_g[:,1] / (gfdl_tclw_lat_g[:,1] + gfdl_tciw_lat_g[:,1])) * gfdl_tcc_lat_g[:,1]

gfdl_tclw_frac_lat_g = np.vstack((gfdl_tclw_lat_g[:,0], gfdl_tclw_frac_lat_g)).T
gfdl_tciw_frac_lat_g = np.vstack((gfdl_tciw_lat_g[:,0], gfdl_tciw_frac_lat_g)).T


#---gfdl Global Profile---#

gfdl_tcc_alt_g = b['cf'][:26] # 0-1
gfdl_tclw_alt_g = b['lw'][:21] #kg/kg
gfdl_tciw_alt_g = b['iw'][:25] #kg/kg
gfdl_temp_alt_g = b['temp'][:] #K
gfdl_plevel_alt_g = b['pressure'][:] #hPa

gfdl_tcc_temp_g = b['cf_t'][:29] # 0-1
gfdl_tclw_temp_g = b['lw_t'][:21] #kg/kg
gfdl_tciw_temp_g = b['iw_t'][:23] #kg/kg

#gfdl_tcc_plevel_g = b['Cloud Fraction with Pressure'][:] # 0-1
#gfdl_tclw_plevel_g = b['Specific Liquid Water Content with Pressure'][:21] #kg/kg
#gfdl_tciw_plevel_g = b['Specific Ice Water Content with Pressure'][:] #kg/kg

#---gfdl Southern Ocean Profile---#

gfdl_tcc_alt_so = b['cf_so'][:23] # 0-1
gfdl_tclw_alt_so = b['lw_so'][:19] #kg/kg
gfdl_tciw_alt_so = b['iw_so'][:23] #kg/kg
gfdl_temp_alt_so = b['temp_so'][:] #K
gfdl_plevel_alt_so = b['pressure_so'][:] #hPa

gfdl_tcc_temp_so = b['cf_t_so'][:23] # 0-1
gfdl_tclw_temp_so = b['lw_t_so'][:19] #kg/kg
gfdl_tciw_temp_so = b['iw_t_so'][:23] #kg/kg

#gfdl_tcc_plevel_so = b['Cloud Fraction with Pressure'][:23] # 0-1
#gfdl_tclw_plevel_so = b['Specific Liquid Water Content with Pressure'][:19] #kg/kg
#gfdl_tciw_plevel_so = b['Specific Ice Water Content with Pressure'][:23] #kg/kg

############################################################################### CAM5 Data

#---CAM5 Global Latitude Data---#

cam_tcc_lat_g = d['tcc'][:] # 0-1
cam_tclw_lat_g = d['tclw'][:] #kgm^-2
cam_tciw_lat_g = d['tciw'][:] #kgm^-2

#---CAM5 Southern Ocean Latitude Data---#

cam_tcc_lat_so = cam_tcc_lat_g[cam_tcc_lat_g[:,0]>=-70]
cam_tcc_lat_so = cam_tcc_lat_so[cam_tcc_lat_so[:,0]<=-50] # 0-1

cam_tclw_lat_so = cam_tclw_lat_g[cam_tclw_lat_g[:,0]>=-70]
cam_tclw_lat_so = cam_tclw_lat_so[cam_tclw_lat_so[:,0]<=-50] #kgm^-2

cam_tciw_lat_so = cam_tciw_lat_g[cam_tciw_lat_g[:,0]>=-70]
cam_tciw_lat_so = cam_tciw_lat_so[cam_tciw_lat_so[:,0]<=-50] #kgm^-2

#---CAM5 Phase Fractions---#

cam_tclw_frac_lat_g = (cam_tclw_lat_g[:,1] / (cam_tclw_lat_g[:,1] + cam_tciw_lat_g[:,1])) * cam_tcc_lat_g[:,1]
cam_tciw_frac_lat_g = (cam_tciw_lat_g[:,1] / (cam_tclw_lat_g[:,1] + cam_tciw_lat_g[:,1])) * cam_tcc_lat_g[:,1]

cam_tclw_frac_lat_g = np.vstack((cam_tclw_lat_g[:,0], cam_tclw_frac_lat_g)).T
cam_tciw_frac_lat_g = np.vstack((cam_tciw_lat_g[:,0], cam_tciw_frac_lat_g)).T

#---CAM5 lat-alt contour data---#

cam_tclw_alt_lat = d['lw_alt_lat'][:] #kg/kg
cam_tciw_alt_lat = d['iw_alt_lat'][:] #kg/kg
cam_lat = d['lat'][:]
cam_alt = d['alt'][:]

#---CAM5 Global Profile---#

#cam_tcc_alt_g = d['cf'][8:] # 0-1
cam_tclw_alt_g = d['lw'][:] #kg/kg
cam_tciw_alt_g = d['iw'][:] #kg/kg
cam_temp_alt_g = d['temp'][:] #K
cam_plevel_alt_g = d['pressure'][:] #hPa

#cam_tcc_temp_g = d['cf_t'][13:] # 0-1
cam_tclw_temp_g = d['lw_t'][:] #kg/kg
cam_tciw_temp_g = d['iw_t'][:] #kg/kg

#cam_tcc_plevel_g = d['Cloud Fraction with Pressure'][8:37] # 0-1
#cam_tclw_plevel_g = d['Specific Liquid Water Content with Pressure'][16:37] #kg/kg
#cam_tciw_plevel_g = d['Specific Ice Water Content with Pressure'][10:37] #kg/kg

#---CAM5 Southern Ocean Profile---#

#cam_tcc_alt_so = d['cf_so'][11:] # 0-1
cam_tclw_alt_so = d['lw_so'][:] #kg/kg
cam_tciw_alt_so = d['iw_so'][:] #kg/kg
cam_temp_alt_so = d['temp_so'][:] #K
#cam_plevel_alt_so = d['pressure_so'][:] #hPa

#cam_tcc_temp_so = d['cf_t_so'][11:] # 0-1
cam_tclw_temp_so = d['lw_t_so'][:] #kg/kg
cam_tciw_temp_so = d['iw_t_so'][:] #kg/kg

#cam_tcc_plevel_so = d['Cloud Fraction with Pressure'][12:37] # 0-1
#cam_tclw_plevel_so = d['Specific Liquid Water Content with Pressure'][18:37] #kg/kg
#cam_tciw_plevel_so = d['Specific Ice Water Content with Pressure'][11:37] #kg/kg

############################################################################### End importing Data

end = time.time()
print('Importing data took:', end - start, 's')

os.chdir('E:/University/University/MSc/Models/Images')


############################################################################### Global Latitude Plots

#---Plot Global Cloud Fraction with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tcc_lat_g[:,0],cccm_tcc_lat_g[:,1], '-r', label='CCCM')
ax.plot(ecmwf_tcc_lat_g[:,0],ecmwf_tcc_lat_g[:,1], '-b', label='ECMWF')
ax.plot(gfdl_tcc_lat_g[:,0],gfdl_tcc_lat_g[:,1], '-g', label='GFDL')
ax.plot(cam_tcc_lat_g[:,0],cam_tcc_lat_g[:,1], '-m', label='CAM5')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 Global Cloud Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tcc_lat_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Global LWP with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_lat_g[:,0],cccm_tclw_lat_g[:,1], '-r', label='CCCM')
ax.plot(ecmwf_tclw_lat_g[:,0],ecmwf_tclw_lat_g[:,1], '-b', label='ECMWF')
ax.plot(gfdl_tclw_lat_g[:,0],gfdl_tclw_lat_g[:,1], '-g', label='GFDL')
ax.plot(cam_tclw_lat_g[:,0],cam_tclw_lat_g[:,1], '-m', label='CAM5')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Content $kgm^{-2}$')

plt.title('07.2006 to 04.2011 Global Liquid Water Content vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_lat_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global IWP with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tciw_lat_g[:,0],cccm_tciw_lat_g[:,1], '--r', label='CCCM')
ax.plot(ecmwf_tciw_lat_g[:,0],ecmwf_tciw_lat_g[:,1], '--b', label='ECMWF')
ax.plot(gfdl_tciw_lat_g[:,0],gfdl_tciw_lat_g[:,1], '--g', label='GFDL')
ax.plot(cam_tciw_lat_g[:,0],cam_tciw_lat_g[:,1], '--m', label='CAM5')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Content $kgm^{-2}$')

plt.title('07.2006 to 04.2011 Global Ice Water Content vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tciw_lat_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Comparison LWP and IWP with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_lat_g[:,0],cccm_tclw_lat_g[:,1], '-r', label='Liquid - CCCM')
ax.plot(ecmwf_tclw_lat_g[:,0],ecmwf_tclw_lat_g[:,1], '-b', label='Liquid - ECMWF')
ax.plot(gfdl_tclw_lat_g[:,0],gfdl_tclw_lat_g[:,1], '-g', label='Liquid - GFDL')
ax.plot(cam_tclw_lat_g[:,0],cam_tclw_lat_g[:,1], '-m', label='Liquid - CAM5')

ax.plot(cccm_tciw_lat_g[:,0],cccm_tciw_lat_g[:,1], '--r', label='Ice - CCCM')
ax.plot(ecmwf_tciw_lat_g[:,0],ecmwf_tciw_lat_g[:,1], '--b', label='Ice - ECMWF')
ax.plot(gfdl_tciw_lat_g[:,0],gfdl_tciw_lat_g[:,1], '--g', label='Ice - GFDL')
ax.plot(cam_tciw_lat_g[:,0],cam_tciw_lat_g[:,1], '--m', label='Ice - CAM5')


ax.legend(loc='upper center', bbox_to_anchor=(1.2, 1.0));

ax.set_xlabel('Latitude')
ax.set_ylabel('LWP and IWP $kgm^{-2}$')

plt.title('07.2006 to 04.2011 Global Liquid and Ice Water Paths vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_tciw_lat_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Global Comparison Ice and Liquid Cloud Fractions with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_frac_lat_g[:,0],cccm_tclw_frac_lat_g[:,1], '-r', label='Liquid - CCCM')
ax.plot(ecmwf_tclw_frac_lat_g[:,0],ecmwf_tclw_frac_lat_g[:,1], '-b', label='Liquid - ECMWF')
ax.plot(gfdl_tclw_frac_lat_g[:,0],gfdl_tclw_frac_lat_g[:,1], '-g', label='Liquid - GFDL')
ax.plot(cam_tclw_frac_lat_g[:,0],cam_tclw_frac_lat_g[:,1], '-m', label='Liquid - CAM5')

ax.plot(cccm_tciw_frac_lat_g[:,0],cccm_tciw_frac_lat_g[:,1], '--r', label='Ice - CCCM')
ax.plot(ecmwf_tciw_frac_lat_g[:,0],ecmwf_tciw_frac_lat_g[:,1], '--b', label='Ice - ECMWF')
ax.plot(gfdl_tciw_frac_lat_g[:,0],gfdl_tciw_frac_lat_g[:,1], '--g', label='Ice - GFDL')
ax.plot(cam_tciw_frac_lat_g[:,0],cam_tciw_frac_lat_g[:,1], '--m', label='Ice - CAM5')


ax.legend(loc='upper center', bbox_to_anchor=(1.2, 1.0));

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Water Fractions')

plt.title('07.2006 to 04.2011 Global Liquid and Ice Water Fractions of Clouds vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_tciw_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CCCM Cloud Fraction and Phase Fraction with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(cccm_tcc_lat_g[:,0],cccm_tcc_lat_g[:,1], '-r', label='Cloud Fraction')
ax1.plot(cccm_tclw_frac_lat_g[:,0],cccm_tclw_frac_lat_g[:,1], '-b', label='Liquid Fraction')
ax1.plot(cccm_tciw_frac_lat_g[:,0],cccm_tciw_frac_lat_g[:,1], '--b', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CCCM Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_CCCM_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot ECMWF Cloud Fraction and Phase Content with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(ecmwf_tcc_lat_g[:,0],ecmwf_tcc_lat_g[:,1], '-r', label='Cloud Fraction')
ax1.plot(ecmwf_tclw_frac_lat_g[:,0],ecmwf_tclw_frac_lat_g[:,1], '-b', label='Liquid Fraction')
ax1.plot(ecmwf_tciw_frac_lat_g[:,0],ecmwf_tciw_frac_lat_g[:,1], '--b', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 ECMWF Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_ECMWF_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""
#---Plot GFDL Cloud Fraction and Phase Fraction with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(gfdl_tcc_lat_g[:,0],gfdl_tcc_lat_g[:,1], '-r', label='Cloud Fraction')
ax1.plot(gfdl_tclw_frac_lat_g[:,0],gfdl_tclw_frac_lat_g[:,1], '-b', label='Liquid Fraction')
ax1.plot(gfdl_tciw_frac_lat_g[:,0],gfdl_tciw_frac_lat_g[:,1], '--b', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 GFDL Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_GFDL_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CAM5 Cloud Fraction and Phase Fraction with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(cam_tcc_lat_g[:,0],cam_tcc_lat_g[:,1], '-r', label='Cloud Fraction')
ax1.plot(cam_tclw_frac_lat_g[:,0],cam_tclw_frac_lat_g[:,1], '-b', label='Liquid Fraction')
ax1.plot(cam_tciw_frac_lat_g[:,0],cam_tciw_frac_lat_g[:,1], '--b', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CAM5 Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_CAM5_lat_frac_g.svg", format="svg", bbox_inches='tight')
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
ax.plot(cam_tcc_lat_so[:,0],cam_tcc_lat_so[:,1], '-m', label='CAM5')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 Southern Ocean Cloud Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tcc_lat_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Southern Ocean Liquid Water Path with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_lat_so[:,0],cccm_tclw_lat_so[:,1], '-r', label='CCCM')
ax.plot(ecmwf_tclw_lat_so[:,0],ecmwf_tclw_lat_so[:,1], '-b', label='ECMWF')
ax.plot(gfdl_tclw_lat_so[:,0],gfdl_tclw_lat_so[:,1], '-g', label='GFDL')
ax.plot(cam_tclw_lat_so[:,0],cam_tclw_lat_so[:,1], '-m', label='CAM5')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Latitude')
ax.set_ylabel('LWP $kgm^{-2}$')

plt.title('07.2006 to 04.2011 Southern Ocean Liquid Water Path vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_lat_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Southern Ocean Ice Water Path with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tciw_lat_so[:,0],cccm_tciw_lat_so[:,1], '--r', label='CCCM')
ax.plot(ecmwf_tciw_lat_so[:,0],ecmwf_tciw_lat_so[:,1], '--b', label='ECMWF')
ax.plot(gfdl_tciw_lat_so[:,0],gfdl_tciw_lat_so[:,1], '--g', label='GFDL')
ax.plot(cam_tciw_lat_so[:,0],cam_tciw_lat_so[:,1], '--m', label='CAM5')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
          ncol=4, fancybox=True, shadow=True);

ax.set_xlabel('Latitude')
ax.set_ylabel('IWP $kgm^{-2}$')

plt.title('07.2006 to 04.2011 Southern Ocean Ice Water Path vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tciw_lat_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Comparison Liquid and Ice Water Path with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_lat_so[:,0],cccm_tclw_lat_so[:,1], '-r', label='Liquid - CCCM')
ax.plot(ecmwf_tclw_lat_so[:,0],ecmwf_tclw_lat_so[:,1], '-b', label='Liquid - ECMWF')
ax.plot(gfdl_tclw_lat_so[:,0],gfdl_tclw_lat_so[:,1], '-g', label='Liquid - GFDL')
ax.plot(cam_tclw_lat_so[:,0],cam_tclw_lat_so[:,1], '-m', label='Liquid - CAM5')
ax.plot(cccm_tciw_lat_so[:,0],cccm_tciw_lat_so[:,1], '--r', label='Ice - CCCM')
ax.plot(ecmwf_tciw_lat_so[:,0],ecmwf_tciw_lat_so[:,1], '--b', label='Ice - ECMWF')
ax.plot(gfdl_tciw_lat_so[:,0],gfdl_tciw_lat_so[:,1], '--g', label='Ice - GFDL')
ax.plot(cam_tciw_lat_so[:,0],cam_tciw_lat_so[:,1], '--m', label='Ice - CAM5')


ax.legend(loc='upper center', bbox_to_anchor=(1.2, 1.0));

ax.set_xlabel('Latitude')
ax.set_ylabel('IWP and LWP $kgm^{-2}$')

plt.title('07.2006 to 04.2011 Southern Ocean Liquid and Ice Water Path vs Latitude')

plt.grid(True)

plt.savefig("07.2006_04.2011_tclw_tciw_lat_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Comparison Ice and Liquid Cloud Fractions with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_frac_lat_so[:,0],cccm_tclw_frac_lat_so[:,1], '-r', label='Liquid - CCCM')
ax.plot(ecmwf_tclw_frac_lat_so[:,0],ecmwf_tclw_frac_lat_so[:,1], '-b', label='Liquid - ECMWF')
ax.plot(gfdl_tclw_frac_lat_so[:,0],gfdl_tclw_frac_lat_so[:,1], '-g', label='Liquid - GFDL')
ax.plot(cam_tclw_frac_lat_so[:,0],cam_tclw_frac_lat_so[:,1], '-m', label='Liquid - CAM5')

ax.plot(cccm_tciw_frac_lat_so[:,0],cccm_tciw_frac_lat_so[:,1], '--r', label='Ice - CCCM')
ax.plot(ecmwf_tciw_frac_lat_so[:,0],ecmwf_tciw_frac_lat_so[:,1], '--b', label='Ice - ECMWF')
ax.plot(gfdl_tciw_frac_lat_so[:,0],gfdl_tciw_frac_lat_so[:,1], '--g', label='Ice - GFDL')
ax.plot(cam_tciw_frac_lat_so[:,0],cam_tciw_frac_lat_so[:,1], '--m', label='Ice - CAM5')


ax.legend(loc='upper center', bbox_to_anchor=(1.2, 1.0));

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Water Fractions')

plt.title('07.2006 to 04.2011 Southern Ocean Liquid and Ice Water Fractions of Clouds vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_tciw_lat_frac_so.svg", format="svg", bbox_inches='tight')
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
ax2.set_ylabel('Cloud Content $kgm^{-2}$')

plt.title('07.2006 to 04.2011 CCCM Southern Ocean Cloud Fraction and Phase Content vs Latitude')


plt.savefig("07.2006_04.2011_CCCM_lat_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot ECMWF Cloud Fraction and Phase Content with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(ecmwf_tcc_lat_so[:,0],ecmwf_tcc_lat_so[:,1], '-r', label='Cloud Fraction')
ax2.plot(ecmwf_tclw_lat_so[:,0],ecmwf_tclw_lat_so[:,1], '-b', label='Liquid Content')
ax2.plot(ecmwf_tciw_lat_so[:,0],ecmwf_tciw_lat_so[:,1], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Cloud Content $kgm^{-2}$')

plt.title('07.2006 to 04.2011 ECMWF Southern Ocean Cloud Fraction and Phase Content vs Latitude')


plt.savefig("07.2006_04.2011_ECMWF_lat_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""
#---Plot GFDL Cloud Fraction and Phase Content with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(gfdl_tcc_lat_so[:,0],gfdl_tcc_lat_so[:,1], '-r', label='Cloud Fraction')
ax2.plot(gfdl_tclw_lat_so[:,0],gfdl_tclw_lat_so[:,1], '-b', label='Liquid Content')
ax2.plot(gfdl_tciw_lat_so[:,0],gfdl_tciw_lat_so[:,1], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Cloud Content $kgm^{-2}$')

plt.title('07.2006 to 04.2011 GFDL Southern Ocean Cloud Fraction and Phase Content vs Latitude')

plt.savefig("07.2006_04.2011_GFDL_lat_so.svg", format="svg", bbox_inches='tight')
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

plt.title('07.2006 to 04.2011 Global Cloud Fraction vs Altitude')

plt.grid(True)

plt.savefig("07.2006_04.2011_tcc_alt_g.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Global Specific Liquid Water Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_alt_g.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Global Specific Ice Water Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tciw_alt_g.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Global Specific Liquid and Ice Water Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_tciw_alt_g.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Global Temperature vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_T_alt_g.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Global Temperature vs Altitude')
#plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("07.2006_04.2011_p_alt_g.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Global Specific Liquid Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_T_g.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Global Specific Ice Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("07.2006_04.2011_tciw_T_g.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Global Specific Liquid and Ice Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_tciw_T_g.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 CCCM Global Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_CCCM_alt_g.svg", format="svg", bbox_inches='tight')
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

plt.title('07.2006 to 04.2011 ECMWF Global Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_ECMWF_alt_g.svg", format="svg", bbox_inches='tight')
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

plt.title('07.2006 to 04.2011 GFDL Global Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_GFDL_alt_g.svg", format="svg", bbox_inches='tight')
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

plt.title('07.2006 to 04.2011 Southern Ocean Cloud Fraction vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tcc_alt_so.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Southern Ocean Specific Liquid Water Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_alt_so.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Southern Ocean Specific Ice Water Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tciw_alt_so.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Southern Ocean Specific Liquid and Ice Water Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_tciw_alt_so.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Southern Ocean Temperature vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_T_alt_so.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Southern Ocean Temperature vs Altitude')
#plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("07.2006_04.2011_P_alt_so.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Southern Ocean Specific Liquid Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_T_so.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Southern Ocean Specific Ice Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("07.2006_04.2011_tciw_T_so.svg", format="svg", bbox_inches='tight')

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

plt.title('07.2006 to 04.2011 Southern Ocean Specific Liquid and Ice Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_tciw_T_so.svg", format="svg", bbox_inches='tight')
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

plt.title('07.2006 to 04.2011 CCCM Southern Ocean Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_CCCM_alt_so.svg", format="svg", bbox_inches='tight')
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

plt.title('07.2006 to 04.2011 ECMWF Southern Ocean Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_ECMWF_alt_so.svg", format="svg", bbox_inches='tight')
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

plt.title('07.2006 to 04.2011 GFDL Southern Ocean Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_GFDL_alt_so.svg", format="svg", bbox_inches='tight')
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

plt.title('07.2006 to 04.2011 Southern Ocean Specific Liquid and Ice Water Content vs Altitude')

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

plt.title('07.2006 to 04.2011 Southern Ocean Specific Liquid and Ice Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.show()
"""