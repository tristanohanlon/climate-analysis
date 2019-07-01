# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Tristan O'Hanlon

This import the reduced datasets from 07.2006 to 04.2011 CCCM, ECMWF and gfdl4. 
The code can select both global or southern ocean data.

"""
import time
import matplotlib.pyplot as plt
import h5py
import os
import numpy as np

start = time.time()


#---Importing Data from Reduced Datasets---#
"""
# Uni Laptop
#ECMWF Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets')
a = h5py.File('07.2006_04.2011_ECMWF.h5', 'r')

#CCCM Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
c = h5py.File('07.2006_04.2011_CCCM.h5', 'r')

#GFDL_AM4 Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GFDL_AM4/reduced_datasets')
b = h5py.File('07.2006_04.2011_GFDL_AM4.h5', 'r')

#GFDL_CM3 Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GFDL_CM3/reduced_datasets')
h = h5py.File('07.2006_04.2011_GFDL_CM3.h5', 'r')

#CAM5 Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CAM5/reduced_datasets')
d = h5py.File('07.2006_04.2011_CAM5.h5', 'r')

#CAM6 Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CAM6/reduced_datasets')
e = h5py.File('07.2006_04.2011_CAM6.h5', 'r')

#CAPLISO-GOCCP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CALIPSO_GOCCP/reduced_datasets')
f = h5py.File('07.2006_04.2011_CALIPSO.h5', 'r')

#CERES Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CERES/reduced_datasets')
g = h5py.File('07.2006_04.2011_CERES.h5', 'r')
"""

# Home PC
#ECMWF Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF/reduced_datasets')
a = h5py.File('07.2006_04.2011_ECMWF.h5', 'r')

#CCCM Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
c = h5py.File('07.2006_04.2011_CCCM.h5', 'r')

#GFDL_AM4 Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL_AM4/reduced_datasets')
b = h5py.File('07.2006_04.2011_GFDL_AM4.h5', 'r')

#GFDL_CM3 Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL_CM3/reduced_datasets')
h = h5py.File('07.2006_04.2011_GFDL_CM3.h5', 'r')

#CAM5 Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CAM5/reduced_datasets')
d = h5py.File('07.2006_04.2011_CAM5.h5', 'r')

#CAM6 Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CAM6/reduced_datasets')
e = h5py.File('07.2006_04.2011_CAM6.h5', 'r')

#CAPLISO-GOCCP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CALIPSO_GOCCP/reduced_datasets')
f = h5py.File('07.2006_04.2011_CALIPSO.h5', 'r')

#CERES Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CERES/reduced_datasets')
g = h5py.File('07.2006_04.2011_CERES.h5', 'r')

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

#CAM6 Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CAM6/reduced_datasets')
e = h5py.File('07.2006_04.2011_CAM6.h5', 'r')

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

ecmwf_tcc_alt_g = a['cf'][:] # 0-1
ecmwf_tclw_alt_g = a['lw'][:] #kg/kg
ecmwf_tciw_alt_g = a['iw'][:] #kg/kg
ecmwf_temp_alt_g = a['temp'][:] #K
ecmwf_plevel_alt_g = a['pressure'][:] #hPa

ecmwf_tcc_temp_g = a['cf_t'][13:] # 0-1
ecmwf_tclw_temp_g = a['lw_t'][17:] #kg/kg
ecmwf_tciw_temp_g = a['iw_t'][10:] #kg/kg

#---ecmwf Phase Profile Fractions---#

ecmwf_tclw_frac_alt_g = (ecmwf_tclw_alt_g[:,1] / (ecmwf_tclw_alt_g[:,1] + ecmwf_tciw_alt_g[:,1])) * ecmwf_tcc_alt_g[:,1]
ecmwf_tciw_frac_alt_g = (ecmwf_tciw_alt_g[:,1] / (ecmwf_tclw_alt_g[:,1] + ecmwf_tciw_alt_g[:,1])) * ecmwf_tcc_alt_g[:,1]

ecmwf_tclw_frac_alt_g = np.vstack((ecmwf_tclw_alt_g[:,0], ecmwf_tclw_frac_alt_g)).T
ecmwf_tciw_frac_alt_g = np.vstack((ecmwf_tciw_alt_g[:,0], ecmwf_tciw_frac_alt_g)).T

ecmwf_tcc_alt_g = a['cf'][8:] # 0-1
ecmwf_tclw_alt_g = a['lw'][14:] #kg/kg
ecmwf_tciw_alt_g = a['iw'][9:] #kg/kg



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

cccm_tcc_alt_g = c['cf'][:] # 0-1 4:101
cccm_tclw_alt_g = c['lw'][:113] #kg/kg 4:63
cccm_tciw_alt_g = c['iw'][:113] #kg/kg 4:93
cccm_temp_alt_g = c['temp'][:] #K
cccm_plevel_alt_g = c['pressure'][:] #hPa

cccm_tcc_temp_g = c['cf_t'][:80] # 0-1
cccm_tclw_temp_g = c['lw_t'][:100] #kg/kg
cccm_tciw_temp_g = c['iw_t'][:100] #kg/kg

#---CCCM Phase Profile Fractions---#

cccm_tclw_frac_alt_g = (cccm_tclw_alt_g[:,1] / (cccm_tclw_alt_g[:,1] + cccm_tciw_alt_g[:,1])) * cccm_tcc_alt_g[:,1]
cccm_tciw_frac_alt_g = (cccm_tciw_alt_g[:,1] / (cccm_tclw_alt_g[:,1] + cccm_tciw_alt_g[:,1])) * cccm_tcc_alt_g[:,1]

cccm_tclw_frac_alt_g = np.vstack((cccm_tclw_alt_g[:,0], cccm_tclw_frac_alt_g)).T
cccm_tciw_frac_alt_g = np.vstack((cccm_tciw_alt_g[:,0], cccm_tciw_frac_alt_g)).T

cccm_tcc_alt_g = c['cf'][4:101] # 0-1
cccm_tclw_alt_g = c['lw'][4:63] #kg/kg
cccm_tciw_alt_g = c['iw'][4:93] #kg/kg



#---CCCM Southern Ocean Profile---#

cccm_tcc_alt_so = c['cf_so'][4:75] # 0-1
cccm_tclw_alt_so = c['lw_so'][4:68] #kg/kg
cccm_tciw_alt_so = c['iw_so'][4:65] #kg/kg
cccm_temp_alt_so = c['temp_so'][:] #K
cccm_plevel_alt_so = c['pressure_so'][:] #hPa

cccm_tcc_temp_so = c['cf_t_so'][:58] # 0-1
cccm_tclw_temp_so = c['lw_t_so'][:100] #kg/kg
cccm_tciw_temp_so = c['iw_t_so'][:100] #kg/kg

#cccm_tcc_plevel_so = c['Cloud Fraction with Pressure'][:] # 0-1
#cccm_tclw_plevel_so = c['Specific Liquid Water Content with Pressure'][:] #kg/kg
#cccm_tciw_plevel_so = c['Specific Ice Water Content with Pressure'][:] #kg/kg

############################################################################### GFDL_AM4 Data

#---GFDL_AM4 Global Latitude Data---#

gfdl4_tcc_lat_g = b['tcc'][:] # 0-1
gfdl4_tclw_lat_g = b['tclw'][:] #kgm^-2
gfdl4_tciw_lat_g = b['tciw'][:] #kgm^-2

#---GFDL_AM4 Southern Ocean Latitude Data---#

gfdl4_tcc_lat_so = gfdl4_tcc_lat_g[gfdl4_tcc_lat_g[:,0]>=-70]
gfdl4_tcc_lat_so = gfdl4_tcc_lat_so[gfdl4_tcc_lat_so[:,0]<=-50] # 0-1

gfdl4_tclw_lat_so = gfdl4_tclw_lat_g[gfdl4_tclw_lat_g[:,0]>=-70]
gfdl4_tclw_lat_so = gfdl4_tclw_lat_so[gfdl4_tclw_lat_so[:,0]<=-50] #kgm^-2

gfdl4_tciw_lat_so = gfdl4_tciw_lat_g[gfdl4_tciw_lat_g[:,0]>=-70]
gfdl4_tciw_lat_so = gfdl4_tciw_lat_so[gfdl4_tciw_lat_so[:,0]<=-50] #kgm^-2

#---GFDL_AM4 Phase Fractions---#

gfdl4_tclw_frac_lat_g = (gfdl4_tclw_lat_g[:,1] / (gfdl4_tclw_lat_g[:,1] + gfdl4_tciw_lat_g[:,1])) * gfdl4_tcc_lat_g[:,1]
gfdl4_tciw_frac_lat_g = (gfdl4_tciw_lat_g[:,1] / (gfdl4_tclw_lat_g[:,1] + gfdl4_tciw_lat_g[:,1])) * gfdl4_tcc_lat_g[:,1]

gfdl4_tclw_frac_lat_g = np.vstack((gfdl4_tclw_lat_g[:,0], gfdl4_tclw_frac_lat_g)).T
gfdl4_tciw_frac_lat_g = np.vstack((gfdl4_tciw_lat_g[:,0], gfdl4_tciw_frac_lat_g)).T

#---GFDL_AM4 lat-alt contour data---#

gfdl4_tclw_alt_lat = b['lw_alt_lat'][:] #kg/kg
gfdl4_tciw_alt_lat = b['iw_alt_lat'][:] #kg/kg
gfdl4_lat = b['lat'][:]
gfdl4_alt = b['alt'][:]


#---GFDL_AM4 Global Profile---#

gfdl4_tcc_alt_g = b['cf'][:] # 0-1
gfdl4_tclw_alt_g = b['lw'][:] #kg/kg
gfdl4_tciw_alt_g = b['iw'][:] #kg/kg
gfdl4_temp_alt_g = b['temp'][:] #K
gfdl4_plevel_alt_g = b['pressure'][:] #hPa

gfdl4_tcc_temp_g = b['cf_t'][:29] # 0-1
gfdl4_tclw_temp_g = b['lw_t'][:21] #kg/kg
gfdl4_tciw_temp_g = b['iw_t'][:23] #kg/kg

#---GFDL_AM4 Phase Profile Fractions---#

gfdl4_tclw_frac_alt_g = (gfdl4_tclw_alt_g[:,1] / (gfdl4_tclw_alt_g[:,1] + gfdl4_tciw_alt_g[:,1])) * gfdl4_tcc_alt_g[:,1]
gfdl4_tciw_frac_alt_g = (gfdl4_tciw_alt_g[:,1] / (gfdl4_tclw_alt_g[:,1] + gfdl4_tciw_alt_g[:,1])) * gfdl4_tcc_alt_g[:,1]

gfdl4_tclw_frac_alt_g = np.vstack((gfdl4_tclw_alt_g[:,0], gfdl4_tclw_frac_alt_g)).T
gfdl4_tciw_frac_alt_g = np.vstack((gfdl4_tciw_alt_g[:,0], gfdl4_tciw_frac_alt_g)).T

gfdl4_tcc_alt_g = b['cf'][:26] # 0-1
gfdl4_tclw_alt_g = b['lw'][:21] #kg/kg
gfdl4_tciw_alt_g = b['iw'][:25] #kg/kg

#gfdl4_tcc_plevel_g = b['Cloud Fraction with Pressure'][:] # 0-1
#gfdl4_tclw_plevel_g = b['Specific Liquid Water Content with Pressure'][:21] #kg/kg
#gfdl4_tciw_plevel_g = b['Specific Ice Water Content with Pressure'][:] #kg/kg

#---GFDL_AM4 Southern Ocean Profile---#

gfdl4_tcc_alt_so = b['cf_so'][:23] # 0-1
gfdl4_tclw_alt_so = b['lw_so'][:19] #kg/kg
gfdl4_tciw_alt_so = b['iw_so'][:23] #kg/kg
gfdl4_temp_alt_so = b['temp_so'][:] #K
gfdl4_plevel_alt_so = b['pressure_so'][:] #hPa

gfdl4_tcc_temp_so = b['cf_t_so'][:23] # 0-1
gfdl4_tclw_temp_so = b['lw_t_so'][:19] #kg/kg
gfdl4_tciw_temp_so = b['iw_t_so'][:23] #kg/kg

#gfdl4_tcc_plevel_so = b['Cloud Fraction with Pressure'][:23] # 0-1
#gfdl4_tclw_plevel_so = b['Specific Liquid Water Content with Pressure'][:19] #kg/kg
#gfdl4_tciw_plevel_so = b['Specific Ice Water Content with Pressure'][:23] #kg/kg


############################################################################### GFDL_CM3 Data

#---GFDL_CM3 Global Latitude Data---#

gfdl3_tcc_lat_g = h['tcc'][:] # 0-1
gfdl3_tclw_lat_g = h['tclw'][:] #kgm^-2
gfdl3_tciw_lat_g = h['tciw'][:] #kgm^-2

#---GFDL_CM3 Southern Ocean Latitude Data---#

gfdl3_tcc_lat_so = gfdl3_tcc_lat_g[gfdl3_tcc_lat_g[:,0]>=-70]
gfdl3_tcc_lat_so = gfdl3_tcc_lat_so[gfdl3_tcc_lat_so[:,0]<=-50] # 0-1

gfdl3_tclw_lat_so = gfdl3_tclw_lat_g[gfdl3_tclw_lat_g[:,0]>=-70]
gfdl3_tclw_lat_so = gfdl3_tclw_lat_so[gfdl3_tclw_lat_so[:,0]<=-50] #kgm^-2

gfdl3_tciw_lat_so = gfdl3_tciw_lat_g[gfdl3_tciw_lat_g[:,0]>=-70]
gfdl3_tciw_lat_so = gfdl3_tciw_lat_so[gfdl3_tciw_lat_so[:,0]<=-50] #kgm^-2

#---GFDL_CM3 Phase Fractions---#

gfdl3_tclw_frac_lat_g = (gfdl3_tclw_lat_g[:,1] / (gfdl3_tclw_lat_g[:,1] + gfdl3_tciw_lat_g[:,1])) * gfdl3_tcc_lat_g[:,1]
gfdl3_tciw_frac_lat_g = (gfdl3_tciw_lat_g[:,1] / (gfdl3_tclw_lat_g[:,1] + gfdl3_tciw_lat_g[:,1])) * gfdl3_tcc_lat_g[:,1]

gfdl3_tclw_frac_lat_g = np.vstack((gfdl3_tclw_lat_g[:,0], gfdl3_tclw_frac_lat_g)).T
gfdl3_tciw_frac_lat_g = np.vstack((gfdl3_tciw_lat_g[:,0], gfdl3_tciw_frac_lat_g)).T

#---GFDL_CM3 lat-alt contour data---#

gfdl3_tclw_alt_lat = h['lw_alt_lat'][:] #kg/kg
gfdl3_tciw_alt_lat = h['iw_alt_lat'][:] #kg/kg
gfdl3_lat = h['lat'][:]
gfdl3_alt = h['alt'][:]


#---GFDL_CM3 Global Profile---#

gfdl3_tcc_alt_g = h['cf'][:] # 0-1
gfdl3_tclw_alt_g = h['lw'][:] #kg/kg
gfdl3_tciw_alt_g = h['iw'][:] #kg/kg
gfdl3_temp_alt_g = h['temp'][:] #K
gfdl3_plevel_alt_g = h['pressure'][:] #hPa

gfdl3_tcc_temp_g = h['cf_t'][:29] # 0-1
gfdl3_tclw_temp_g = h['lw_t'][:21] #kg/kg
gfdl3_tciw_temp_g = h['iw_t'][:23] #kg/kg

#---GFDL_CM3 Phase Profile Fractions---#

gfdl3_tclw_frac_alt_g = (gfdl3_tclw_alt_g[:,1] / (gfdl3_tclw_alt_g[:,1] + gfdl3_tciw_alt_g[:,1])) * gfdl3_tcc_alt_g[:,1]
gfdl3_tciw_frac_alt_g = (gfdl3_tciw_alt_g[:,1] / (gfdl3_tclw_alt_g[:,1] + gfdl3_tciw_alt_g[:,1])) * gfdl3_tcc_alt_g[:,1]

gfdl3_tclw_frac_alt_g = np.vstack((gfdl3_tclw_alt_g[:,0], gfdl3_tclw_frac_alt_g)).T
gfdl3_tciw_frac_alt_g = np.vstack((gfdl3_tciw_alt_g[:,0], gfdl3_tciw_frac_alt_g)).T

gfdl3_tcc_alt_g = b['cf'][:26] # 0-1
gfdl3_tclw_alt_g = b['lw'][:21] #kg/kg
gfdl3_tciw_alt_g = b['iw'][:25] #kg/kg


#gfdl3_tcc_plevel_g = h['Cloud Fraction with Pressure'][:] # 0-1
#gfdl3_tclw_plevel_g = h['Specific Liquid Water Content with Pressure'][:21] #kg/kg
#gfdl3_tciw_plevel_g = h['Specific Ice Water Content with Pressure'][:] #kg/kg

#---GFDL_CM3 Southern Ocean Profile---#

gfdl3_tcc_alt_so = h['cf_so'][:23] # 0-1
gfdl3_tclw_alt_so = h['lw_so'][:19] #kg/kg
gfdl3_tciw_alt_so = h['iw_so'][:23] #kg/kg
gfdl3_temp_alt_so = h['temp_so'][:] #K

gfdl3_tcc_temp_so = h['cf_t_so'][:23] # 0-1
gfdl3_tclw_temp_so = h['lw_t_so'][:19] #kg/kg
gfdl3_tciw_temp_so = h['iw_t_so'][:23] #kg/kg

#gfdl3_tcc_plevel_so = h['Cloud Fraction with Pressure'][:23] # 0-1
#gfdl3_tclw_plevel_so = h['Specific Liquid Water Content with Pressure'][:19] #kg/kg
#gfdl3_tciw_plevel_so = h['Specific Ice Water Content with Pressure'][:23] #kg/kg

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
cam_lat = np.hstack(cam_lat)
cam_alt = d['alt'][:]
cam_alt = np.hstack(cam_alt)


#---CAM5 Global Profile---#

cam_tcc_alt_g = d['cf'][:] # 0-1
cam_tclw_alt_g = d['lw'][:] #kg/kg
cam_tciw_alt_g = d['iw'][:] #kg/kg
cam_temp_alt_g = d['temp'][:] #K
cam_plevel_alt_g = d['pressure'][:] #hPa

cam_tcc_temp_g = d['cf_t'][13:] # 0-1
cam_tclw_temp_g = d['lw_t'][18:29] #kg/kg
cam_tciw_temp_g = d['iw_t'][10:29] #kg/kg

#---CAM5 Phase Profile Fractions---#

cam_tclw_frac_alt_g = (cam_tclw_alt_g[:,1] / (cam_tclw_alt_g[:,1] + cam_tciw_alt_g[:,1])) * cam_tcc_alt_g[:,1]
cam_tciw_frac_alt_g = (cam_tciw_alt_g[:,1] / (cam_tclw_alt_g[:,1] + cam_tciw_alt_g[:,1])) * cam_tcc_alt_g[:,1]

cam_tclw_frac_alt_g = np.vstack((cam_tclw_alt_g[:,0], cam_tclw_frac_alt_g)).T
cam_tciw_frac_alt_g = np.vstack((cam_tciw_alt_g[:,0], cam_tciw_frac_alt_g)).T

cam_tcc_alt_g = d['cf'][8:] # 0-1
cam_tclw_alt_g = d['lw'][13:29] #kg/kg
cam_tciw_alt_g = d['iw'][7:29] #kg/kg

#cam_tcc_plevel_g = d['Cloud Fraction with Pressure'][8:37] # 0-1
#cam_tclw_plevel_g = d['Specific Liquid Water Content with Pressure'][16:37] #kg/kg
#cam_tciw_plevel_g = d['Specific Ice Water Content with Pressure'][10:37] #kg/kg

#---CAM5 Southern Ocean Profile---#

cam_tcc_alt_so = d['cf_so'][11:] # 0-1
cam_tclw_alt_so = d['lw_so'][13:29] #kg/kg
cam_tciw_alt_so = d['iw_so'][10:29] #kg/kg
cam_temp_alt_so = d['temp_so'][19:29] #K
#cam_plevel_alt_so = d['pressure_so'][:] #hPa

cam_tcc_temp_so = d['cf_t_so'][11:] # 0-1
cam_tclw_temp_so = d['lw_t_so'][20:29] #kg/kg
cam_tciw_temp_so = d['iw_t_so'][11:29] #kg/kg

#cam_tcc_plevel_so = d['Cloud Fraction with Pressure'][12:37] # 0-1
#cam_tclw_plevel_so = d['Specific Liquid Water Content with Pressure'][18:37] #kg/kg
#cam_tciw_plevel_so = d['Specific Ice Water Content with Pressure'][11:37] #kg/kg

############################################################################### CAM6 Data

#---CAM6 Global Latitude Data---#

cam6_tcc_lat_g = e['tcc'][:] # 0-1
cam6_tclw_lat_g = e['tclw'][:] #kgm^-2
cam6_tciw_lat_g = e['tciw'][:] #kgm^-2

#---CAM6 Southern Ocean Latitude Data---#

cam6_tcc_lat_so = cam6_tcc_lat_g[cam6_tcc_lat_g[:,0]>=-70]
cam6_tcc_lat_so = cam6_tcc_lat_so[cam6_tcc_lat_so[:,0]<=-50] # 0-1

cam6_tclw_lat_so = cam6_tclw_lat_g[cam6_tclw_lat_g[:,0]>=-70]
cam6_tclw_lat_so = cam6_tclw_lat_so[cam6_tclw_lat_so[:,0]<=-50] #kgm^-2

cam6_tciw_lat_so = cam6_tciw_lat_g[cam6_tciw_lat_g[:,0]>=-70]
cam6_tciw_lat_so = cam6_tciw_lat_so[cam6_tciw_lat_so[:,0]<=-50] #kgm^-2

#---CAM6 Phase Fractions---#

cam6_tclw_frac_lat_g = (cam6_tclw_lat_g[:,1] / (cam6_tclw_lat_g[:,1] + cam6_tciw_lat_g[:,1])) * cam6_tcc_lat_g[:,1]
cam6_tciw_frac_lat_g = (cam6_tciw_lat_g[:,1] / (cam6_tclw_lat_g[:,1] + cam6_tciw_lat_g[:,1])) * cam6_tcc_lat_g[:,1]

cam6_tclw_frac_lat_g = np.vstack((cam6_tclw_lat_g[:,0], cam6_tclw_frac_lat_g)).T
cam6_tciw_frac_lat_g = np.vstack((cam6_tciw_lat_g[:,0], cam6_tciw_frac_lat_g)).T

#---CAM6 lat-alt contour data---#

cam6_tcc_alt_lat = e['cf_alt_lat'][:] 
cam6_tclw_alt_lat = e['lw_alt_lat'][:] #kg/kg
cam6_tciw_alt_lat = e['iw_alt_lat'][:] #kg/kg
cam6_lat = e['lat'][:]
cam6_lat = np.hstack(cam6_lat)
cam6_alt = e['alt'][:]
cam6_alt = np.hstack(cam6_alt)


#---CAM6 Global Profile---#

cam6_tcc_alt_g = e['cf'][:] # 0-1
cam6_tclw_alt_g = e['lw'][:] #kg/kg
cam6_tciw_alt_g = e['iw'][:] #kg/kg
cam6_temp_alt_g = e['temp'][:] #K
cam6_plevel_alt_g = e['pressure'][:] #hPa

cam6_tcc_temp_g = e['cf_t'][13:] # 0-1
cam6_tclw_temp_g = e['lw_t'][18:] #kg/kg
cam6_tciw_temp_g = e['iw_t'][10:] #kg/kg

#---CAM6 Phase Profile Fractions---#

cam6_tclw_frac_alt_g = (cam6_tclw_alt_g[:,1] / (cam6_tclw_alt_g[:,1] + cam6_tciw_alt_g[:,1])) * cam6_tcc_alt_g[:,1]
cam6_tciw_frac_alt_g = (cam6_tciw_alt_g[:,1] / (cam6_tclw_alt_g[:,1] + cam6_tciw_alt_g[:,1])) * cam6_tcc_alt_g[:,1]

cam6_tclw_frac_alt_g = np.vstack((cam6_tclw_alt_g[:,0], cam6_tclw_frac_alt_g)).T
cam6_tciw_frac_alt_g = np.vstack((cam6_tciw_alt_g[:,0], cam6_tciw_frac_alt_g)).T

cam6_tcc_alt_g = e['cf'][8:] # 0-1
cam6_tclw_alt_g = e['lw'][13:] #kg/kg
cam6_tciw_alt_g = e['iw'][7:] #kg/kg


#cam6_tcc_plevel_g = e['Cloud Fraction with Pressure'][8:37] # 0-1
#cam6_tclw_plevel_g = e['Specific Liquid Water Content with Pressure'][16:37] #kg/kg
#cam6_tciw_plevel_g = e['Specific Ice Water Content with Pressure'][10:37] #kg/kg

#---CAM6 Southern Ocean Profile---#

cam6_tcc_alt_so = e['cf_so'][11:] # 0-1
cam6_tclw_alt_so = e['lw_so'][13:] #kg/kg
cam6_tciw_alt_so = e['iw_so'][10:] #kg/kg
cam6_temp_alt_so = e['temp_so'][19:] #K
#cam6_plevel_alt_so = e['pressure_so'][:] #hPa

cam6_tcc_temp_so = e['cf_t_so'][11:] # 0-1
cam6_tclw_temp_so = e['lw_t_so'][20:] #kg/kg
cam6_tciw_temp_so = e['iw_t_so'][11:] #kg/kg

#cam6_tcc_plevel_so = e['Cloud Fraction with Pressure'][12:37] # 0-1
#cam6_tclw_plevel_so = e['Specific Liquid Water Content with Pressure'][18:37] #kg/kg
#cam6_tciw_plevel_so = e['Specific Ice Water Content with Pressure'][11:37] #kg/kg


############################################################################### CALIPSO_GOCCP Data

#---CALIPSO_GOCCP Global Latitude Data---#

calipso_tcc_lat_g = f['tcc'][:] # 0-1
calipso_tclw_frac_lat_g = f['tclw_frac'][:] # 0-1
calipso_tciw_frac_lat_g = f['tciw_frac'][:] # 0-1

#---CALIPSO_GOCCP Southern Ocean Latitude Data---#

calipso_tcc_lat_so = calipso_tcc_lat_g[calipso_tcc_lat_g[:,0]>=-70]
calipso_tcc_lat_so = calipso_tcc_lat_so[calipso_tcc_lat_so[:,0]<=-50] # 0-1

calipso_tclw_frac_lat_so = calipso_tclw_frac_lat_g[calipso_tclw_frac_lat_g[:,0]>=-70]
calipso_tclw_frac_lat_so = calipso_tclw_frac_lat_so[calipso_tclw_frac_lat_so[:,0]<=-50] # 0-1

calipso_tciw_frac_lat_so = calipso_tciw_frac_lat_g[calipso_tciw_frac_lat_g[:,0]>=-70]
calipso_tciw_frac_lat_so = calipso_tciw_frac_lat_so[calipso_tciw_frac_lat_so[:,0]<=-50] # 0-1


#---CALIPSO_GOCCP Global Profile---#

calipso_tcc_alt_g = f['cf'][:] # 0-1
calipso_tclw_frac_alt_g = f['lw_frac'][:] # 0-1
calipso_tciw_frac_alt_g = f['iw_frac'][:] # 0-1

############################################################################### CERES Data

#---CERES Global Latitude Data---#

ceres_tcc_lat_g = g['tcc'][:] # 0-1
ceres_tclw_lat_g = g['tclw'][:] #kgm^-2
ceres_tciw_lat_g = g['tciw'][:] #kgm^-2
ceres_tclw_frac_lat_g = g['tclw_frac'][:] # 0-1
ceres_tciw_frac_lat_g = g['tciw_frac'][:] # 0-1

#---CERES Southern Ocean Latitude Data---#

ceres_tcc_lat_so = ceres_tcc_lat_g[ceres_tcc_lat_g[:,0]>=-70]
ceres_tcc_lat_so = ceres_tcc_lat_so[ceres_tcc_lat_so[:,0]<=-50] # 0-1

ceres_tclw_lat_so = ceres_tclw_lat_g[ceres_tclw_lat_g[:,0]>=-70]
ceres_tclw_lat_so = ceres_tclw_lat_so[ceres_tclw_lat_so[:,0]<=-50] #kgm^-2

ceres_tciw_lat_so = ceres_tciw_lat_g[ceres_tciw_lat_g[:,0]>=-70]
ceres_tciw_lat_so = ceres_tciw_lat_so[ceres_tciw_lat_so[:,0]<=-50] #kgm^-2


############################################################################### End importing Data

end = time.time()
print('Importing data took:', end - start, 's')

#os.chdir('E:/University/University/MSc/Models/Images')
os.chdir('E:/University/University/MSc/Models/Images/Meeting 20.6')

############################################################################### Contour plots

#---Plot GFDL_AM4 Specific Liquid Water Content---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.contourf(gfdl4_lat, gfdl4_alt[0:18], gfdl4_tclw_alt_lat[0:18])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')

plt.title('07.2006 to 04.2011 GFDL_AM4 Specific Liquid Water Content')

plt.savefig("07.2006_04.2011_contour_gfdl4_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot GFDL_AM4 Specific Ice Water Content---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.contourf(gfdl4_lat, gfdl4_alt[:26], gfdl4_tciw_alt_lat[:26])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')

plt.title('07.2006 to 04.2011 GFDL_AM4 Specific Ice Water Content')

plt.savefig("07.2006_04.2011_contour_gfdl4_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""
#---Plot GFDL_CM3 Specific Liquid Water Content---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.contourf(gfdl3_lat, gfdl3_alt[0:18], gfdl3_tclw_alt_lat[0:18])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')

plt.title('07.2006 to 04.2011 GFDL_CM3 Specific Liquid Water Content')

plt.savefig("07.2006_04.2011_contour_gfdl3_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot GFDL_CM3 Specific Ice Water Content---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.contourf(gfdl3_lat, gfdl3_alt[:26], gfdl3_tciw_alt_lat[:26])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')

plt.title('07.2006 to 04.2011 GFDL_CM3 Specific Ice Water Content')

plt.savefig("07.2006_04.2011_contour_gfdl3_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""
#---Plot CAM5 Specific Liquid Water Content---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.contourf(cam_lat, cam_alt[20:30], cam_tclw_alt_lat[20:30])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')

plt.title('07.2006 to 04.2011 CAM5 Specific Liquid Water Content')

plt.savefig("07.2006_04.2011_contour_CAM5_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CAM5 Specific Ice Water Content---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.contourf(cam_lat, cam_alt[8:30], cam_tciw_alt_lat[8:30])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')

plt.title('07.2006 to 04.2011 CAM5 Specific Ice Water Content')

plt.savefig("07.2006_04.2011_contour_CAM5_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CAM6 Specific Liquid Water Content---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.contourf(cam6_lat, cam6_alt[22:32], cam6_tclw_alt_lat[22:32])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')

plt.title('07.2006 to 04.2011 CAM6 Specific Liquid Water Content')

plt.savefig("07.2006_04.2011_contour_CAM6_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CAM6 Specific Ice Water Content---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.contourf(cam6_lat, cam6_alt[10:32], cam6_tciw_alt_lat[10:32])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')

plt.title('07.2006 to 04.2011 CAM6 Specific Ice Water Content')

plt.savefig("07.2006_04.2011_contour_CAM6_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

############################################################################### Global Latitude Plots

#---Plot Global Cloud Fraction with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

#ax.plot(cccm_tcc_lat_g[:,0],cccm_tcc_lat_g[:,1], '-r', label='CCCM')
ax.plot(calipso_tcc_lat_g[:,0],calipso_tcc_lat_g[:,1], '-r', label='CAPLISO')
ax.plot(ceres_tcc_lat_g[:,0],ceres_tcc_lat_g[:,1], '-b', label='CERES')

ax.plot(ecmwf_tcc_lat_g[:,0],ecmwf_tcc_lat_g[:,1], '-k', label='ECMWF-ERA5')
ax.plot(gfdl4_tcc_lat_g[:,0],gfdl4_tcc_lat_g[:,1], '-g', label='CMIP6-GFDL-AM4')
ax.plot(gfdl3_tcc_lat_g[:,0],gfdl3_tcc_lat_g[:,1], '-y', label='CMIP5-GFDL-CM3')
ax.plot(cam_tcc_lat_g[:,0],cam_tcc_lat_g[:,1], '-m', label='CMIP5-CESM1-CAM5-RPC4.5')
ax.plot(cam6_tcc_lat_g[:,0],cam6_tcc_lat_g[:,1], '-c', label='CMIP6-CESM2.1-CAM6-RPC4.5')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

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

#ax.plot(cccm_tclw_lat_g[:,0],cccm_tclw_lat_g[:,1], ':r', label='CCCM')
#ax.plot(calipso_tclw_lat_g[:,0],calipso_tclw_lat_g[:,1], ':', color='brown', label='CAPLISO-GOCCP')
#ax.plot(ceres_tclw_lat_g[:,0],ceres_tclw_lat_g[:,1], ':b', label='CERES')

#ax.plot(ecmwf_tclw_lat_g[:,0],ecmwf_tclw_lat_g[:,1], ':b', label='ECMWF')
#ax.plot(gfdl4_tclw_lat_g[:,0],gfdl4_tclw_lat_g[:,1], ':g', label='CMIP6-GFDL-AM4')
#ax.plot(gfdl3_tclw_lat_g[:,0],gfdl3_tclw_lat_g[:,1], ':y', label='CMIP5-GFDL-CM3')
#ax.plot(cam_tclw_lat_g[:,0],cam_tclw_lat_g[:,1], ':m', label='CMIP5-CESM1-CAM5-RPC4.5')
#ax.plot(cam6_tclw_lat_g[:,0],cam6_tclw_lat_g[:,1], ':c', label='CMIP6-CESM2.1-CAM6-RPC4.5')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=5);

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
#ax.plot(calipso_tciw_lat_g[:,0],calipso_tciw_lat_g[:,1], '--', color='orange', label='CAPLISO-GOCCP')
ax.plot(ceres_tciw_lat_g[:,0],ceres_tciw_lat_g[:,1], '--b', label='CERES')

#ax.plot(ecmwf_tciw_lat_g[:,0],ecmwf_tciw_lat_g[:,1], '--b', label='ECMWF')
#ax.plot(gfdl4_tciw_lat_g[:,0],gfdl4_tciw_lat_g[:,1], '--g', label='CMIP6-GFDL-AM4')
#ax.plot(gfdl3_tciw_lat_g[:,0],gfdl3_tciw_lat_g[:,1], '-y', label='CMIP5-GFDL-CM3')
#ax.plot(cam_tciw_lat_g[:,0],cam_tciw_lat_g[:,1], '--m', label='CMIP5-CESM1-CAM5-RPC4.5')
#ax.plot(cam6_tciw_lat_g[:,0],cam6_tciw_lat_g[:,1], '--c', label='CMIP6-CESM2.1-CAM6')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=5);

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

ax.plot(cccm_tclw_lat_g[:,0],cccm_tclw_lat_g[:,1], ':r', label='Liquid - CCCM')
ax.plot(ecmwf_tclw_lat_g[:,0],ecmwf_tclw_lat_g[:,1], ':b', label='Liquid - ECMWF')
ax.plot(gfdl4_tclw_lat_g[:,0],gfdl4_tclw_lat_g[:,1], ':g', label='Liquid - CMIP6-GFDL-AM4')
ax.plot(cam_tclw_lat_g[:,0],cam_tclw_lat_g[:,1], ':m', label='Liquid - CMIP5-CESM1-CAM5-RPC4.5')
ax.plot(cam6_tclw_lat_g[:,0],cam6_tclw_lat_g[:,1], ':c', label='CMIP6-CESM2.1-CAM6')

ax.plot(cccm_tciw_lat_g[:,0],cccm_tciw_lat_g[:,1], '--r', label='Ice - CCCM')
ax.plot(ecmwf_tciw_lat_g[:,0],ecmwf_tciw_lat_g[:,1], '--b', label='Ice - ECMWF')
ax.plot(gfdl4_tciw_lat_g[:,0],gfdl4_tciw_lat_g[:,1], '--g', label='Ice - CMIP6-GFDL-AM4')
ax.plot(cam_tciw_lat_g[:,0],cam_tciw_lat_g[:,1], '--m', label='Ice - CMIP5-CESM1-CAM5-RPC4.5')
ax.plot(cam6_tciw_lat_g[:,0],cam6_tciw_lat_g[:,1], '--c', label='CMIP6-CESM2.1-CAM6')


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

ax.plot(calipso_tclw_frac_lat_g[:,0],calipso_tclw_frac_lat_g[:,1], ':r', label='Liquid - CALIPSO')
ax.plot(ecmwf_tclw_frac_lat_g[:,0],ecmwf_tclw_frac_lat_g[:,1], ':k', label='Liquid - ECMWF-ERA5')
ax.plot(gfdl3_tclw_frac_lat_g[:,0],gfdl3_tclw_frac_lat_g[:,1], ':b', label='Liquid - CMIP5-GFDL-CM3')
ax.plot(gfdl4_tclw_frac_lat_g[:,0],gfdl4_tclw_frac_lat_g[:,1], ':g', label='Liquid - CMIP6-GFDL-AM4')
ax.plot(cam_tclw_frac_lat_g[:,0],cam_tclw_frac_lat_g[:,1], ':m', label='Liquid - CMIP5-CESM1-CAM5-RPC4.5')
ax.plot(cam6_tclw_frac_lat_g[:,0],cam6_tclw_frac_lat_g[:,1], ':c', label='Liquid - CMIP6-CESM2.1-CAM6')

#ax.plot(calipso_tciw_frac_lat_g[:,0],calipso_tciw_frac_lat_g[:,1], '--r', label='Ice - CALIPSO')
#ax.plot(ecmwf_tciw_frac_lat_g[:,0],ecmwf_tciw_frac_lat_g[:,1], '--k', label='Ice - ECMWF-ERA5')
#ax.plot(gfdl3_tciw_frac_lat_g[:,0],gfdl3_tciw_frac_lat_g[:,1], '--b', label='Ice - CMIP5-GFDL-CM3')
#ax.plot(gfdl4_tciw_frac_lat_g[:,0],gfdl4_tciw_frac_lat_g[:,1], '--g', label='Ice - CMIP6-GFDL-AM4')
#ax.plot(cam_tciw_frac_lat_g[:,0],cam_tciw_frac_lat_g[:,1], '--m', label='Ice - CMIP5-CESM1-CAM5-RPC4.5')
#ax.plot(cam6_tciw_frac_lat_g[:,0],cam6_tciw_frac_lat_g[:,1], '--c', label='Ice - CMIP6-CESM2.1-CAM6')


ax.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Liquid Water Fraction')

plt.title('07.2006 to 04.2011 Global Liquid Fractions of Clouds vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CCCM Cloud Fraction and Phase Fraction with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(cccm_tcc_lat_g[:,0],cccm_tcc_lat_g[:,1], '-r', label='Cloud Fraction')
ax1.plot(cccm_tclw_frac_lat_g[:,0],cccm_tclw_frac_lat_g[:,1], ':r', label='Liquid Fraction')
ax1.plot(cccm_tciw_frac_lat_g[:,0],cccm_tciw_frac_lat_g[:,1], '--r', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CCCM Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_CCCM_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CALIPSO-GOCCP Cloud Fraction and Phase Fraction with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(calipso_tcc_lat_g[:,0],calipso_tcc_lat_g[:,1], '-b', label='Cloud Fraction')
ax1.plot(calipso_tclw_frac_lat_g[:,0],calipso_tclw_frac_lat_g[:,1], ':b', label='Liquid Fraction')
ax1.plot(calipso_tciw_frac_lat_g[:,0],calipso_tciw_frac_lat_g[:,1], '--b', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CALIPSO-GOCCP Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_CALIPSO-GOCCP_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CERES Cloud Fraction and Phase Fraction with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(ceres_tcc_lat_g[:,0],ceres_tcc_lat_g[:,1], '-b', label='Cloud Fraction')
ax1.plot(ceres_tclw_frac_lat_g[:,0],ceres_tclw_frac_lat_g[:,1], ':b', label='Liquid Fraction')
ax1.plot(ceres_tciw_frac_lat_g[:,0],ceres_tciw_frac_lat_g[:,1], '--b', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CERES Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_CERES_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Combined Satellite Cloud Fraction and Phase Fraction with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()
ax1.plot(cccm_tcc_lat_g[:,0],cccm_tcc_lat_g[:,1], '-r', label='CCCM Cloud Fraction')
ax1.plot(cccm_tclw_frac_lat_g[:,0],cccm_tclw_frac_lat_g[:,1], ':r', label='CCCM Liquid Fraction')
ax1.plot(cccm_tciw_frac_lat_g[:,0],cccm_tciw_frac_lat_g[:,1], '--r', label='CCCM Ice Fraction')
ax1.plot(calipso_tcc_lat_g[:,0],calipso_tcc_lat_g[:,1], '-', color='orange', label='CALIPSO-GOCCP Cloud Fraction')
ax1.plot(calipso_tclw_frac_lat_g[:,0],calipso_tclw_frac_lat_g[:,1], ':', color='orange', label='CALIPSO-GOCCP Liquid Fraction')
ax1.plot(calipso_tciw_frac_lat_g[:,0],calipso_tciw_frac_lat_g[:,1], '--', color='orange', label='CALIPSO-GOCCP Ice Fraction')
ax1.plot(ceres_tcc_lat_g[:,0],ceres_tcc_lat_g[:,1], '-b', label='CERES Cloud Fraction')
ax1.plot(ceres_tclw_frac_lat_g[:,0],ceres_tclw_frac_lat_g[:,1], ':b', label='CERES Liquid Fraction')
ax1.plot(ceres_tciw_frac_lat_g[:,0],ceres_tciw_frac_lat_g[:,1], '--b', label='CERES Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 Satellite Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_SAT_lat_frac_g.svg", format="svg", bbox_inches='tight')
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
#---Plot CMIP6-GFDL-AM4 Cloud Fraction and Phase Fraction with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(gfdl4_tcc_lat_g[:,0],gfdl4_tcc_lat_g[:,1], '-r', label='Cloud Fraction')
ax1.plot(gfdl4_tclw_frac_lat_g[:,0],gfdl4_tclw_frac_lat_g[:,1], '-b', label='Liquid Fraction')
ax1.plot(gfdl4_tciw_frac_lat_g[:,0],gfdl4_tciw_frac_lat_g[:,1], '--b', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CMIP6-GFDL-AM4 Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_gfdl4_lat_frac_g.svg", format="svg", bbox_inches='tight')
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

plt.title('07.2006 to 04.2011 CMIP5-CESM1-CAM5-RPC4.5 Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_CAM5_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CAM6 Cloud Fraction and Phase Fraction with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(cam6_tcc_lat_g[:,0],cam6_tcc_lat_g[:,1], '-r', label='Cloud Fraction')
ax1.plot(cam6_tclw_frac_lat_g[:,0],cam6_tclw_frac_lat_g[:,1], '-b', label='Liquid Fraction')
ax1.plot(cam6_tciw_frac_lat_g[:,0],cam6_tciw_frac_lat_g[:,1], '--b', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CMIP6-CESM2.1-CAM6 Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_CAM6_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CAM5 and CAM6 Cloud Fraction and Phase Fraction with Latitude Comparisons---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(cam_tcc_lat_g[:,0],cam_tcc_lat_g[:,1], '-r', label='CAM5-RPC4.5 - Cloud Fraction')
ax1.plot(cam_tclw_frac_lat_g[:,0],cam_tclw_frac_lat_g[:,1], '-b', label='CAM5-RPC4.5 - Liquid Fraction')
ax1.plot(cam_tciw_frac_lat_g[:,0],cam_tciw_frac_lat_g[:,1], '--b', label='CAM5-RPC4.5 - Ice Fraction')

ax1.plot(cam6_tcc_lat_g[:,0],cam6_tcc_lat_g[:,1], '-m', label='CAM6 - Cloud Fraction')
ax1.plot(cam6_tclw_frac_lat_g[:,0],cam6_tclw_frac_lat_g[:,1], '-g', label='CAM6 - Liquid Fraction')
ax1.plot(cam6_tciw_frac_lat_g[:,0],cam6_tciw_frac_lat_g[:,1], '--g', label='CAM6 - Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CMIP6-CESM2.1-CAM6 and CMIP5-CESM1-CAM5-RPC4.5 Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_CAM5_CAM6_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""
############################################################################### Southern Ocean Latitude Plots

#---Plot Southern Ocean Cloud Fraction with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tcc_lat_so[:,0],cccm_tcc_lat_so[:,1], '-r', label='CCCM')
ax.plot(ecmwf_tcc_lat_so[:,0],ecmwf_tcc_lat_so[:,1], '-b', label='ECMWF')
ax.plot(gfdl4_tcc_lat_so[:,0],gfdl4_tcc_lat_so[:,1], '-g', label='gfdl4')
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
ax.plot(gfdl4_tclw_lat_so[:,0],gfdl4_tclw_lat_so[:,1], '-g', label='gfdl4')
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
ax.plot(gfdl4_tciw_lat_so[:,0],gfdl4_tciw_lat_so[:,1], '--g', label='gfdl4')
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
ax.plot(gfdl4_tclw_lat_so[:,0],gfdl4_tclw_lat_so[:,1], '-g', label='Liquid - gfdl4')
ax.plot(cam_tclw_lat_so[:,0],cam_tclw_lat_so[:,1], '-m', label='Liquid - CAM5')
ax.plot(cccm_tciw_lat_so[:,0],cccm_tciw_lat_so[:,1], '--r', label='Ice - CCCM')
ax.plot(ecmwf_tciw_lat_so[:,0],ecmwf_tciw_lat_so[:,1], '--b', label='Ice - ECMWF')
ax.plot(gfdl4_tciw_lat_so[:,0],gfdl4_tciw_lat_so[:,1], '--g', label='Ice - gfdl4')
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
ax.plot(gfdl4_tclw_frac_lat_so[:,0],gfdl4_tclw_frac_lat_so[:,1], '-g', label='Liquid - gfdl4')
ax.plot(cam_tclw_frac_lat_so[:,0],cam_tclw_frac_lat_so[:,1], '-m', label='Liquid - CAM5')

ax.plot(cccm_tciw_frac_lat_so[:,0],cccm_tciw_frac_lat_so[:,1], '--r', label='Ice - CCCM')
ax.plot(ecmwf_tciw_frac_lat_so[:,0],ecmwf_tciw_frac_lat_so[:,1], '--b', label='Ice - ECMWF')
ax.plot(gfdl4_tciw_frac_lat_so[:,0],gfdl4_tciw_frac_lat_so[:,1], '--g', label='Ice - gfdl4')
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
#---Plot gfdl4 Cloud Fraction and Phase Content with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(gfdl4_tcc_lat_so[:,0],gfdl4_tcc_lat_so[:,1], '-r', label='Cloud Fraction')
ax2.plot(gfdl4_tclw_lat_so[:,0],gfdl4_tclw_lat_so[:,1], '-b', label='Liquid Content')
ax2.plot(gfdl4_tciw_lat_so[:,0],gfdl4_tciw_lat_so[:,1], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Cloud Content $kgm^{-2}$')

plt.title('07.2006 to 04.2011 gfdl4 Southern Ocean Cloud Fraction and Phase Content vs Latitude')

plt.savefig("07.2006_04.2011_gfdl4_lat_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CAM5 Cloud Fraction and Phase Content with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(cam_tcc_lat_so[:,0],cam_tcc_lat_so[:,1], '-r', label='Cloud Fraction')
ax2.plot(cam_tclw_lat_so[:,0],cam_tclw_lat_so[:,1], '-b', label='Liquid Content')
ax2.plot(cam_tciw_lat_so[:,0],cam_tciw_lat_so[:,1], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')
ax2.set_ylabel('Cloud Content $kgm^{-2}$')

plt.title('07.2006 to 04.2011 CAM5 Southern Ocean Cloud Fraction and Phase Content vs Latitude')

plt.savefig("07.2006_04.2011_CAM5_lat_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""



############################################################################### Global Altitude Plots

#---Plot Global Cloud Fraction Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

#ax.plot(cccm_tcc_alt_g[:,1],cccm_tcc_alt_g[:,0], '-r', label='CCCM')
ax.plot(calipso_tcc_alt_g[:,1],calipso_tcc_alt_g[:,0], '-r', label='CALIPSO')

ax.plot(ecmwf_tcc_alt_g[:,1],ecmwf_tcc_alt_g[:,0], '-k', label='ECMWF-ERA5')
ax.plot(gfdl3_tcc_alt_g[:,1],gfdl3_tcc_alt_g[:,0], '-b', label='CMIP5-GFDL-CM3')
ax.plot(gfdl4_tcc_alt_g[:,1],gfdl4_tcc_alt_g[:,0], '-g', label='CMIP6-GFDL-AM4')
ax.plot(cam_tcc_alt_g[:,1],cam_tcc_alt_g[:,0], '-m', label='CMIP5-CESM1-CAM5-RPC4.5')
ax.plot(cam6_tcc_alt_g[:,1],cam6_tcc_alt_g[:,0], '-c', label='CMIP6-CESM2.1-CAM6-RPC4.5')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

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

ax.plot(calipso_tclw_alt_g[:,1],calipso_tclw_alt_g[:,0], '-r', label='CALIPSO')
ax.plot(ecmwf_tclw_alt_g[:,1],ecmwf_tclw_alt_g[:,0], '-k', label='ECMWF-ERA5')
ax.plot(gfdl3_tclw_alt_g[:,1],gfdl3_tclw_alt_g[:,0], '-b', label='CMIP5-GFDL-CM3')
ax.plot(gfdl4_tclw_alt_g[:,1],gfdl4_tclw_alt_g[:,0], '-g', label='CMIP6-GFDL-AM4')
ax.plot(cam_tclw_alt_g[:,1],cam_tclw_alt_g[:,0], '-m', label='CMIP5-CESM1-CAM5-RPC4.5')
ax.plot(cam6_tclw_alt_g[:,1],cam6_tclw_alt_g[:,0], '-c', label='CMIP6-CESM2.1-CAM6-RPC4.5')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('07.2006 to 04.2011 Global Specific Liquid Water Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""
#---Plot Global Cloud Liquid Water Fraction Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(calipso_tclw_frac_alt_g[:21,1],calipso_tclw_frac_alt_g[:21,0], '-r', label='CALIPSO')
ax.plot(ecmwf_tclw_frac_alt_g[17:36,1],ecmwf_tclw_frac_alt_g[17:36,0], '-k', label='ECMWF-ERA5')
ax.plot(gfdl3_tclw_frac_alt_g[:20,1],gfdl3_tclw_frac_alt_g[:20,0], '-b', label='CMIP5-GFDL-CM3')
ax.plot(gfdl4_tclw_frac_alt_g[:20,1],gfdl4_tclw_frac_alt_g[:20,0], '-g', label='CMIP6-GFDL-AM4')
ax.plot(cam_tclw_frac_alt_g[16:,1],cam_tclw_frac_alt_g[16:,0], '-m', label='CMIP5-CESM1-CAM5-RPC4.5')
ax.plot(cam6_tclw_frac_alt_g[16:,1],cam6_tclw_frac_alt_g[16:,0], '-c', label='CMIP6-CESM2.1-CAM6-RPC4.5')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Cloud Liquid Water Fraction')

plt.title('07.2006 to 04.2011 Cloud Liquid Water Fraction vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_frac_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""
#---Plot Global Specific Ice Water Content Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tciw_alt_g[:,1]*10000,cccm_tciw_alt_g[:,0], '--r', label='CCCM')
ax.plot(ecmwf_tciw_alt_g[:,1]*10000,ecmwf_tciw_alt_g[:,0], '--b', label='ECMWF')
ax.plot(gfdl4_tciw_alt_g[:,1]*10000,gfdl4_tciw_alt_g[:,0], '--g', label='gfdl4')
ax.plot(cam_tciw_alt_g[:,1]*10000,cam_tciw_alt_g[:,0], '--m', label='CAM5')
ax.plot(cam6_tciw_alt_g[:,1]*10000,cam6_tciw_alt_g[:,0], '--c', label='CAM6')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('07.2006 to 04.2011 Global Specific Ice Water Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tciw_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""
#---Plot Global Cloud Ice Water Fraction Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(calipso_tciw_frac_alt_g[:,1],calipso_tciw_frac_alt_g[:,0], '-r', label='CALIPSO')
ax.plot(ecmwf_tciw_frac_alt_g[:,1],ecmwf_tciw_frac_alt_g[:,0], '-k', label='ECMWF-ERA5')
ax.plot(gfdl3_tciw_frac_alt_g[:,1],gfdl3_tciw_frac_alt_g[:,0], '-b', label='CMIP5-GFDL-CM3')
ax.plot(gfdl4_tciw_frac_alt_g[:,1],gfdl4_tciw_frac_alt_g[:,0], '-g', label='CMIP6-GFDL-AM4')
ax.plot(cam_tciw_frac_alt_g[:,1],cam_tciw_frac_alt_g[:,0], '-m', label='CMIP5-CESM1-CAM5-RPC4.5')
ax.plot(cam6_tciw_frac_alt_g[:,1],cam6_tciw_frac_alt_g[:,0], '-c', label='CMIP6-CESM2.1-CAM6-RPC4.5')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Cloud Ice Water Fraction')

plt.title('07.2006 to 04.2011 Cloud Ice Water Fraction vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tciw_frac_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""
#---Plot Global Comparison Specific Liquid and Ice Water Content Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_alt_g[:,1]*10000,cccm_tclw_alt_g[:,0], '-r', label='Liquid - CCCM')
ax.plot(ecmwf_tclw_alt_g[:,1]*10000,ecmwf_tclw_alt_g[:,0], '-b', label='Liquid - ECMWF')
ax.plot(gfdl4_tclw_alt_g[:,1]*10000,gfdl4_tclw_alt_g[:,0], '-g', label='Liquid - gfdl4')
ax.plot(cam_tclw_alt_g[:,1]*10000,cam_tclw_alt_g[:,0], '-m', label='Liquid - CAM5')
ax.plot(cccm_tciw_alt_g[:,1]*10000,cccm_tciw_alt_g[:,0], '--r', label='Ice - CCCM')
ax.plot(ecmwf_tciw_alt_g[:,1]*10000,ecmwf_tciw_alt_g[:,0], '--b', label='Ice - ECMWF')
ax.plot(gfdl4_tciw_alt_g[:,1]*10000,gfdl4_tciw_alt_g[:,0], '--g', label='Ice - gfdl4')
ax.plot(cam_tciw_alt_g[:,1]*10000,cam_tciw_alt_g[:,0], '--m', label='Ice - CAM5')


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
ax.plot(gfdl4_temp_alt_g[:,1],gfdl4_temp_alt_g[:,0], '-g', label='gfdl4')
ax.plot(cam_temp_alt_g[:,1],cam_temp_alt_g[:,0], '-m', label='CAM5')
ax.plot(cam6_temp_alt_g[:,1],cam6_temp_alt_g[:,0], '-c', label='CAM6')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

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
ax.plot(gfdl4_plevel_alt_g[:,1],gfdl4_plevel_alt_g[:,0], '-g', label='gfdl4')
ax.plot(cam_plevel_alt_g[:,1],cam_plevel_alt_g[:,0], '-m', label='CAM5')

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
ax.plot(gfdl4_tclw_temp_g[:,1]*10000,gfdl4_tclw_temp_g[:,0], '-g', label='gfdl4')
ax.plot(cam_tclw_temp_g[:,1]*10000,cam_tclw_temp_g[:,0], '-m', label='CAM5')
ax.plot(cam6_tclw_temp_g[:,1]*10000,cam6_tclw_temp_g[:,0], '-c', label='CAM5')

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
ax.plot(gfdl4_tciw_temp_g[:,1]*10000,gfdl4_tciw_temp_g[:,0], '--g', label='gfdl4')
ax.plot(cam_tciw_temp_g[:,1]*10000,cam_tciw_temp_g[:,0], '--m', label='CAM5')
ax.plot(cam6_tciw_temp_g[:,1]*10000,cam6_tciw_temp_g[:,0], '--c', label='CAM6')

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

ax.plot(cccm_tclw_temp_g[:,1]*10000,cccm_tclw_temp_g[:,0], ':r', label='Liquid - CCCM')
#ax.plot(ecmwf_tclw_temp_g[:,1]*10000,ecmwf_tclw_temp_g[:,0], '-b', label='Liquid - ECMWF')
#ax.plot(gfdl4_tclw_temp_g[:,1]*10000,gfdl4_tclw_temp_g[:,0], '-g', label='Liquid - gfdl4')
#ax.plot(cam_tclw_temp_g[:,1]*10000,cam_tclw_temp_g[:,0], '-m', label='Liquid - CAM5')

ax.plot(cccm_tciw_temp_g[:,1]*10000,cccm_tciw_temp_g[:,0], '--r', label='Ice - CCCM')
#ax.plot(ecmwf_tciw_temp_g[:,1]*10000,ecmwf_tciw_temp_g[:,0], '--b', label='Ice - ECMWF')
#ax.plot(gfdl4_tciw_temp_g[:,1]*10000,gfdl4_tciw_temp_g[:,0], '--g', label='Ice - gfdl4')
#ax.plot(cam_tciw_temp_g[:,1]*10000,cam_tciw_temp_g[:,0], '--m', label='Ice - CAM5')


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
ax2.plot(cccm_tclw_alt_g[:,1]*10000,cccm_tclw_alt_g[:,0], ':r', label='Liquid Content')
ax2.plot(cccm_tciw_alt_g[:,1]*10000,cccm_tciw_alt_g[:,0], '--r', label='Ice Content')

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

#---Plot CALIPSO-GOCCP Cloud Fraction and Phase Profile---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(calipso_tcc_alt_g[:,1],calipso_tcc_alt_g[:,0], '-b', label='Cloud Fraction')
ax1.plot(calipso_tclw_frac_alt_g[:,1],calipso_tclw_frac_alt_g[:,0], ':b', label='Liquid Content')
ax1.plot(calipso_tciw_frac_alt_g[:,1],calipso_tciw_frac_alt_g[:,0], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_ylabel('Altitude')
ax1.set_xlabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CCCM Global Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_CALIPSO-GOCCP_alt_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""


#---Plot Satellite Cloud Fraction and Phase Profile---#

"""
plt.figure()
fig, ax1 = plt.subplots()
ax1.plot(cccm_tcc_alt_g[:,1],cccm_tcc_alt_g[:,0], '-r', label='CCCM Cloud Fraction')
ax1.plot(cccm_tclw_frac_alt_g[:,1],cccm_tclw_frac_alt_g[:,0], ':r', label='CCCM Liquid Content')
ax1.plot(cccm_tciw_frac_alt_g[:,1],cccm_tciw_frac_alt_g[:,0], '--r', label='CCCM Ice Content')

ax1.plot(calipso_tcc_alt_g[:,1],calipso_tcc_alt_g[:,0], '-b', label='CALIPSO-GOCCP Cloud Fraction')
ax1.plot(calipso_tclw_frac_alt_g[:,1],calipso_tclw_frac_alt_g[:,0], ':b', label='CALIPSO-GOCCP Liquid Content')
ax1.plot(calipso_tciw_frac_alt_g[:,1],calipso_tciw_frac_alt_g[:,0], '--b', label='CALIPSO-GOCCP Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax1.set_ylabel('Altitude')
ax1.set_xlabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 Satellite Global Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_SAT_alt_g.svg", format="svg", bbox_inches='tight')
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

#---Plot gfdl4 Cloud Fraction and Phase Profile---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()
ax1.plot(gfdl4_tcc_alt_g[:,1],gfdl4_tcc_alt_g[:,0], '-r', label='Cloud Fraction')
ax2.plot(gfdl4_tclw_alt_g[:,1]*10000,gfdl4_tclw_alt_g[:,0], '-b', label='Liquid Content')
ax2.plot(gfdl4_tciw_alt_g[:,1]*10000,gfdl4_tciw_alt_g[:,0], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_ylabel('Altitude')
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('07.2006 to 04.2011 gfdl4 Global Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_gfdl4_alt_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CAM5 Cloud Fraction and Phase Profile---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(cam_tclw_alt_g[:,1]*10000,cam_tclw_alt_g[:,0], '-b', label='Liquid Content')
ax1.plot(cam_tciw_alt_g[:,1]*10000,cam_tciw_alt_g[:,0], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_ylabel('Altitude')
ax1.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('07.2006 to 04.2011 CAM5 Global Cloud Phase Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_CAM5_alt_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""
############################################################################### Southern Ocean Altitude Plots

#---Plot Southern Ocean Cloud Fraction Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tcc_alt_so[:,1],cccm_tcc_alt_so[:,0], '-r', label='CCCM')
ax.plot(ecmwf_tcc_alt_so[:,1],ecmwf_tcc_alt_so[:,0], '-b', label='ECMWF')
ax.plot(gfdl4_tcc_alt_so[:,1],gfdl4_tcc_alt_so[:,0], '-g', label='gfdl4')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

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
ax.plot(gfdl4_tclw_alt_so[:,1]*10000,gfdl4_tclw_alt_so[:,0], '-g', label='gfdl4')
ax.plot(cam_tclw_alt_so[:,1]*10000,cam_tclw_alt_so[:,0], '-m', label='CAM5')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

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
ax.plot(gfdl4_tciw_alt_so[:,1]*10000,gfdl4_tciw_alt_so[:,0], '--g', label='gfdl4')
ax.plot(cam_tciw_alt_so[:,1]*10000,cam_tciw_alt_so[:,0], '--m', label='CAM5')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

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
ax.plot(gfdl4_tclw_alt_so[:,1]*10000,gfdl4_tclw_alt_so[:,0], '-g', label='Liquid - gfdl4')
ax.plot(cam_tclw_alt_so[:,1]*10000,cam_tclw_alt_so[:,0], '-m', label='Liquid - CAM5')

ax.plot(cccm_tciw_alt_so[:,1]*10000,cccm_tciw_alt_so[:,0], '--r', label='Ice - CCCM')
ax.plot(ecmwf_tciw_alt_so[:,1]*10000,ecmwf_tciw_alt_so[:,0], '--b', label='Ice - ECMWF')
ax.plot(gfdl4_tciw_alt_so[:,1]*10000,gfdl4_tciw_alt_so[:,0], '--g', label='Ice - gfdl4')
ax.plot(cam_tciw_alt_so[:,1]*10000,cam_tciw_alt_so[:,0], '--m', label='Ice - CAM5')


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
ax.plot(gfdl4_temp_alt_so[:,1],gfdl4_temp_alt_so[:,0], '-g', label='gfdl4')
ax.plot(cam_temp_alt_so[:,1],cam_temp_alt_so[:,0], '-m', label='CAM5')

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
ax.plot(gfdl4_plevel_alt_so[:,1],gfdl4_plevel_alt_so[:,0], '-g', label='gfdl4')

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
ax.plot(gfdl4_tclw_temp_so[:,1]*10000,gfdl4_tclw_temp_so[:,0], '-g', label='gfdl4')
ax.plot(cam_tclw_temp_so[:,1]*10000,cam_tclw_temp_so[:,0], '-m', label='CAM5')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

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
ax.plot(gfdl4_tciw_temp_so[:,1]*10000,gfdl4_tciw_temp_so[:,0], '--g', label='gfdl4')
ax.plot(cam_tciw_temp_so[:,1]*10000,cam_tciw_temp_so[:,0], '--m', label='CAM5')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

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


ax.plot(cccm_tclw_temp_so[:,1]*10000,cccm_tclw_temp_so[:,0], ':r', label='Liquid - CCCM')
#ax.plot(ecmwf_tclw_temp_so[:,1]*10000,ecmwf_tclw_temp_so[:,0], '-b', label='Liquid - ECMWF')
#ax.plot(gfdl4_tclw_temp_so[:,1]*10000,gfdl4_tclw_temp_so[:,0], '-g', label='Liquid - gfdl4')
#ax.plot(cam_tclw_temp_so[:,1]*10000,cam_tclw_temp_so[:,0], '-m', label='Liquid - CAM5')

ax.plot(cccm_tciw_temp_so[:,1]*10000,cccm_tciw_temp_so[:,0], '--r', label='Ice - CCCM')
#ax.plot(ecmwf_tciw_temp_so[:,1]*10000,ecmwf_tciw_temp_so[:,0], '--b', label='Ice - ECMWF')
#ax.plot(gfdl4_tciw_temp_so[:,1]*10000,gfdl4_tciw_temp_so[:,0], '--g', label='Ice - gfdl4')
#ax.plot(cam_tciw_temp_so[:,1]*10000,cam_tciw_temp_so[:,0], '--m', label='Ice - CAM5')

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

#---Plot gfdl4 Cloud Fraction and Phase Profile---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax2 = ax1.twiny()
ax1.plot(gfdl4_tcc_alt_so[:,1],gfdl4_tcc_alt_so[:,0], '-r', label='Cloud Fraction')
ax2.plot(gfdl4_tclw_alt_so[:,1]*10000,gfdl4_tclw_alt_so[:,0], '-b', label='Liquid Content')
ax2.plot(gfdl4_tciw_alt_so[:,1]*10000,gfdl4_tciw_alt_so[:,0], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));
ax2.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15));

ax1.set_ylabel('Altitude')
ax1.set_xlabel('Cloud Fraction')
ax2.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('07.2006 to 04.2011 gfdl4 Southern Ocean Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_gfdl4_alt_so.svg", format="svg", bbox_inches='tight')
plt.show()
"""


#---Plot Global and Southern Ocean Comparison Specific Liquid and Ice Water Content vs Altitude---#

"""
plt.figure()
fig, ax = plt.subplots()

#ax.plot(cccm_tclw_alt_g[:,1],cccm_tclw_alt_g[:,0], '-r', label='Liquid - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tclw_alt_g[:,1],ecmwf_tclw_alt_g[:,0], '-b', label='Liquid - ECMWF.ERA5 Reanalysis Model')
#ax.plot(gfdl4_tclw_alt_g[:,1],gfdl4_tclw_alt_g[:,0], '-g', label='Liquid - gfdl4')
#ax.plot(cccm_tciw_alt_g[:,1],cccm_tciw_alt_g[:,0], '--r', label='Ice - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tciw_alt_g[:,1],ecmwf_tciw_alt_g[:,0], '--b', label='Ice - ECMWF.ERA5 Reanalysis Model')
#ax.plot(gfdl4_tciw_alt_g[:,1],gfdl4_tciw_alt_g[:,0], '--g', label='Ice - gfdl4')
#ax.plot(cccm_tclw_alt_so[:,1],cccm_tclw_alt_so[:,0], '-r', label='Liquid - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tclw_alt_so[:,1],ecmwf_tclw_alt_so[:,0], '-b', label='Liquid - ECMWF.ERA5 Reanalysis Model')
#ax.plot(gfdl4_tclw_alt_so[:,1],gfdl4_tclw_alt_so[:,0], '-g', label='Liquid - gfdl4')
#ax.plot(cccm_tciw_alt_so[:,1],cccm_tciw_alt_so[:,0], '--r', label='Ice - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tciw_alt_so[:,1],ecmwf_tciw_alt_so[:,0], '--b', label='Ice - ECMWF.ERA5 Reanalysis Model')
#ax.plot(gfdl4_tciw_alt_so[:,1],gfdl4_tciw_alt_so[:,0], '--g', label='Ice - gfdl4')

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
#ax.plot(gfdl4_tclw_temp_g[:,1],gfdl4_tclw_temp_g[:,0], '-g', label='Liquid - gfdl4')
#ax.plot(cccm_tciw_temp_g[:,1],cccm_tciw_temp_g[:,0], '--r', label='Ice - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tciw_temp_g[:,1],ecmwf_tciw_temp_g[:,0], '--b', label='Ice - ECMWF.ERA5 Reanalysis Model')
#ax.plot(gfdl4_tciw_temp_g[:,1],gfdl4_tciw_temp_g[:,0], '--g', label='Ice - gfdl4')
#ax.plot(cccm_tclw_temp_so[:,1],cccm_tclw_temp_so[:,0], '-r', label='Liquid - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tclw_temp_so[:,1],ecmwf_tclw_temp_so[:,0], '-b', label='Liquid - ECMWF.ERA5 Reanalysis Model')
#ax.plot(gfdl4_tclw_temp_so[:,1],gfdl4_tclw_temp_so[:,0], '-g', label='Liquid - gfdl4')
#ax.plot(cccm_tciw_temp_so[:,1],cccm_tciw_temp_so[:,0], '--r', label='Ice - CCCM Merged Satellite Dataset')
ax.plot(ecmwf_tciw_temp_so[:,1],ecmwf_tciw_temp_so[:,0], '--b', label='Ice - ECMWF.ERA5 Reanalysis Model')
#ax.plot(gfdl4_tciw_temp_so[:,1],gfdl4_tciw_temp_so[:,0], '--g', label='Ice - gfdl4')

ax.legend(loc='upper center', bbox_to_anchor=(1.2, 1.0));

ax.set_ylabel('Temperature (K)')
ax.set_xlabel('Specific Liquid and Ice Water Content $(kg/kg) x 10^{-4}$')

plt.title('07.2006 to 04.2011 Southern Ocean Specific Liquid and Ice Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.show()
"""