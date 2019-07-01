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
import scipy.interpolate as interp

start = time.time()


#---Importing Data from Reduced Datasets---#
"""
# Uni Laptop
#ECMWF-ERA5 Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/ECMWF-ERA5/reduced_datasets')
a = h5py.File('07.2006_04.2011_ecmwf.h5', 'r')

#GFDL-AM4-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/GFDL-AM4-AMIP/reduced_datasets')
b = h5py.File('07.2006_04.2011_gfdl_am4.h5', 'r')

#CCCM Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
c = h5py.File('07.2006_04.2011_cccm.h5', 'r')

#MRI-ESM2-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/MRI-ESM2-AMIP/reduced_datasets')
d = h5py.File('07.2006_04.2011_mri_esm2.h5', 'r')

#CESM2-CAM6-AMIP Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CESM2-CAM6-AMIP/reduced_datasets')
e = h5py.File('07.2006_04.2011_cesm2_cam6.h5', 'r')

#CAPLISO Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CALIPSO_GOCCP/reduced_datasets')
f = h5py.File('07.2006_04.2011_calipso.h5', 'r')

#CERES Data
os.chdir('C:/Users/toha006/University/University/MSc/Models/climate-analysis/CERES/reduced_datasets')
g = h5py.File('07.2006_04.2011_ceres.h5', 'r')
"""

# Home PC
#ECMWF Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/ECMWF-ERA5/reduced_datasets')
a = h5py.File('07.2006_04.2011_ecmwf.h5', 'r')

#GFDL-AM4-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/GFDL-AM4-AMIP/reduced_datasets')
b = h5py.File('07.2006_04.2011_gfdl_am4.h5', 'r')

#CCCM Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
c = h5py.File('07.2006_04.2011_cccm.h5', 'r')

#MRI-ESM2-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/MRI-ESM2-AMIP/reduced_datasets')
d = h5py.File('07.2006_04.2011_mri_esm2.h5', 'r')

#CESM2-CAM6-AMIP Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CESM2-CAM6-AMIP/reduced_datasets')
e = h5py.File('07.2006_04.2011_cesm2_cam6.h5', 'r')

#CAPLISO Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CALIPSO_GOCCP/reduced_datasets')
f = h5py.File('07.2006_04.2011_calipso.h5', 'r')

#CERES Data
os.chdir('E:/University/University/MSc/Models/climate-analysis/CERES/reduced_datasets')
g = h5py.File('07.2006_04.2011_ceres.h5', 'r')

"""
# Laptop
#ECMWF Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/ECMWF-ERA5/reduced_datasets')
a = h5py.File('07.2006_04.2011_ecmwf.h5', 'r')

#GFDL-AM4-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/GFDL-AM4-AMIP/reduced_datasets')
b = h5py.File('07.2006_04.2011_gfdl_am4.h5', 'r')

#CCCM Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CCCM/reduced_datasets')
c = h5py.File('07.2006_04.2011_cccm.h5', 'r')

#MRI-ESM2-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/MRI-ESM2-AMIP/reduced_datasets')
d = h5py.File('07.2006_04.2011_mri_esm2.h5', 'r')

#CESM2-CAM6-AMIP Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CESM2-CAM6-AMIP/reduced_datasets')
e = h5py.File('07.2006_04.2011_cesm2_cam6.h5', 'r')

#CAPLISO Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CALIPSO_GOCCP/reduced_datasets')
f = h5py.File('07.2006_04.2011_calipso.h5', 'r')

#CERES Data
os.chdir('C:/Users/tristan/University/University/MSc/Models/climate-analysis/CERES/reduced_datasets')
g = h5py.File('07.2006_04.2011_ceres.h5', 'r')
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

#---ECMWF lat-alt contour data---#

ecmwf_tcc_alt_lat = b['cf_alt_lat'][:] #kg/kg
ecmwf_tclw_alt_lat = b['lw_alt_lat'][:] #kg/kg
ecmwf_tciw_alt_lat = b['iw_alt_lat'][:] #kg/kg
ecmwf_temp_alt_lat = b['temp_alt_lat'][:] #kg/kg
ecmwf_lat = b['lat'][:]
ecmwf_alt = b['alt'][:]

ecmwf_tclw_frac_alt_lat = (ecmwf_tclw_alt_lat / (ecmwf_tclw_alt_lat + ecmwf_tciw_alt_lat)) * (ecmwf_tcc_alt_lat)
ecmwf_tciw_frac_alt_lat = (ecmwf_tciw_alt_lat / (ecmwf_tclw_alt_lat + ecmwf_tciw_alt_lat)) * (ecmwf_tcc_alt_lat)


#---ECMWF Global Profile---#

ecmwf_tcc_alt_g = a['cf'][:] # 0-1
ecmwf_tclw_alt_g = a['lw'][:] #kg/kg
ecmwf_tciw_alt_g = a['iw'][:] #kg/kg
ecmwf_temp_alt_g = a['temp'][:] #K
ecmwf_plevel_alt_g = a['pressure'][:] #hPa

ecmwf_tcc_temp_g = a['cf_t'][:] # 0-1
ecmwf_tclw_temp_g = a['lw_t'][:] #kg/kg
ecmwf_tciw_temp_g = a['iw_t'][:] #kg/kg

ecmwf_tclw_frac_temp_g = (ecmwf_tclw_temp_g[:,1] / (ecmwf_tclw_temp_g[:,1] + ecmwf_tciw_temp_g[:,1])) * ecmwf_tcc_temp_g[:,1]
ecmwf_tciw_frac_temp_g = (ecmwf_tciw_temp_g[:,1] / (ecmwf_tclw_temp_g[:,1] + ecmwf_tciw_temp_g[:,1])) * ecmwf_tcc_temp_g[:,1]

ecmwf_tclw_frac_temp_g = np.vstack((ecmwf_tclw_temp_g[:,0], ecmwf_tclw_frac_temp_g)).T
ecmwf_tciw_frac_temp_g = np.vstack((ecmwf_tciw_temp_g[:,0], ecmwf_tciw_frac_temp_g)).T


#---ECMWF Phase Profile Fractions---#

ecmwf_tclw_frac_alt_g = (ecmwf_tclw_alt_g[:,1] / (ecmwf_tclw_alt_g[:,1] + ecmwf_tciw_alt_g[:,1])) * ecmwf_tcc_alt_g[:,1]
ecmwf_tciw_frac_alt_g = (ecmwf_tciw_alt_g[:,1] / (ecmwf_tclw_alt_g[:,1] + ecmwf_tciw_alt_g[:,1])) * ecmwf_tcc_alt_g[:,1]

ecmwf_tclw_frac_alt_g = np.vstack((ecmwf_tclw_alt_g[:,0], ecmwf_tclw_frac_alt_g)).T
ecmwf_tciw_frac_alt_g = np.vstack((ecmwf_tciw_alt_g[:,0], ecmwf_tciw_frac_alt_g)).T

ecmwf_tcc_alt_g = a['cf'][8:] # 0-1
ecmwf_tclw_alt_g = a['lw'][14:] #kg/kg
ecmwf_tciw_alt_g = a['iw'][9:] #kg/kg


#---ECMWF Southern Ocean Profile---#

ecmwf_tcc_alt_so = a['cf_so'][11:] # 0-1
ecmwf_tclw_alt_so = a['lw_so'][18:] #kg/kg
ecmwf_tciw_alt_so = a['iw_so'][11:] #kg/kg
ecmwf_temp_alt_so = a['temp_so'][:] #K
ecmwf_plevel_alt_so = a['pressure_so'][:] #hPa

ecmwf_tcc_temp_so = a['cf_t_so'][:] # 0-1
ecmwf_tclw_temp_so = a['lw_t_so'][:] #kg/kg
ecmwf_tciw_temp_so = a['iw_t_so'][:] #kg/kg


ecmwf_tclw_frac_temp_so = (ecmwf_tclw_temp_so[:,1] / (ecmwf_tclw_temp_so[:,1] + ecmwf_tciw_temp_so[:,1])) * ecmwf_tcc_temp_so[:,1]
ecmwf_tciw_frac_temp_so = (ecmwf_tciw_temp_so[:,1] / (ecmwf_tclw_temp_so[:,1] + ecmwf_tciw_temp_so[:,1])) * ecmwf_tcc_temp_so[:,1]

ecmwf_tclw_frac_temp_so = np.vstack((ecmwf_tclw_temp_so[:,0], ecmwf_tclw_frac_temp_so)).T
ecmwf_tciw_frac_temp_so = np.vstack((ecmwf_tciw_temp_so[:,0], ecmwf_tciw_frac_temp_so)).T


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

cccm_tclw_frac_lat_g = c['tclw_frac'][:]
cccm_tciw_frac_lat_g = c['tciw_frac'][:]

#---CCCM Global Profile---#

cccm_tcc_alt_g = c['cf'][:] # 0-1 4:101
cccm_tclw_alt_g = c['lw'][:113] #kg/kg 4:63
cccm_tciw_alt_g = c['iw'][:113] #kg/kg 4:93
cccm_temp_alt_g = c['temp'][:] #K
cccm_plevel_alt_g = c['pressure'][:] #hPa

cccm_tcc_temp_g = c['cf_t'][:] # 0-1
cccm_tclw_temp_g = c['lw_t'][:] #kg/kg
cccm_tciw_temp_g = c['iw_t'][:] #kg/kg

cccm_tclw_frac_alt_g = c['lw_frac'][:]
cccm_tciw_frac_alt_g = c['iw_frac'][:]

#---CCCM Temperature Phase Fractions---#

cccm_tclw_frac_temp_g = c['lw_frac_temp'][:]
cccm_tciw_frac_temp_g = c['iw_frac_temp'][:]

cccm_tclw_frac_temp_so = c['lw_frac_temp_so'][:]
cccm_tciw_frac_temp_so = c['iw_frac_temp_so'][:]


#---CCCM Southern Ocean Profile---#

cccm_tcc_alt_so = c['cf_so'][:] # 0-1
cccm_tclw_alt_so = c['lw_so'][:] #kg/kg
cccm_tciw_alt_so = c['iw_so'][:] #kg/kg
cccm_temp_alt_so = c['temp_so'][:] #K
cccm_plevel_alt_so = c['pressure_so'][:] #hPa

cccm_tcc_temp_so = c['cf_t_so'][:] # 0-1
cccm_tclw_temp_so = c['lw_t_so'][:] #kg/kg
cccm_tciw_temp_so = c['iw_t_so'][:] #kg/kg

cccm_tclw_frac_alt_so = c['lw_frac_so'][:]
cccm_tciw_frac_alt_so = c['iw_frac_so'][:]





############################################################################### GFDL-AM4-AMIP Data

#---GFDL-AM4-AMIP Global Latitude Data---#

gfdl4_tcc_lat_g = b['tcc'][:] # 0-1
gfdl4_tclw_lat_g = b['tclw'][:] #kgm^-2
gfdl4_tciw_lat_g = b['tciw'][:] #kgm^-2

#---GFDL-AM4-AMIP Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
gfdl4_tclw_frac_lat_g = b['tclw_frac'][:] 
gfdl4_tciw_frac_lat_g = b['tciw_frac'][:]


#---GFDL-AM4-AMIP Southern Ocean Latitude Data---#

gfdl4_tcc_lat_so = gfdl4_tcc_lat_g[gfdl4_tcc_lat_g[:,0]>=-70]
gfdl4_tcc_lat_so = gfdl4_tcc_lat_so[gfdl4_tcc_lat_so[:,0]<=-50] # 0-1

gfdl4_tclw_lat_so = gfdl4_tclw_lat_g[gfdl4_tclw_lat_g[:,0]>=-70]
gfdl4_tclw_lat_so = gfdl4_tclw_lat_so[gfdl4_tclw_lat_so[:,0]<=-50] #kgm^-2

gfdl4_tciw_lat_so = gfdl4_tciw_lat_g[gfdl4_tciw_lat_g[:,0]>=-70]
gfdl4_tciw_lat_so = gfdl4_tciw_lat_so[gfdl4_tciw_lat_so[:,0]<=-50] #kgm^-2


#---GFDL-AM4-AMIP lat-alt contour data---#

gfdl4_tcc_alt_lat = b['cf_alt_lat'][:] #kg/kg
gfdl4_tclw_alt_lat = b['lw_alt_lat'][:] #kg/kg
gfdl4_tciw_alt_lat = b['iw_alt_lat'][:] #kg/kg
gfdl4_temp_alt_lat = b['temp_alt_lat'][:] #kg/kg
gfdl4_lat = b['lat'][:]
gfdl4_alt = b['alt'][:]

gfdl4_tclw_frac_alt_lat = (gfdl4_tclw_alt_lat / (gfdl4_tclw_alt_lat + gfdl4_tciw_alt_lat)) * (gfdl4_tcc_alt_lat)
gfdl4_tciw_frac_alt_lat = (gfdl4_tciw_alt_lat / (gfdl4_tclw_alt_lat + gfdl4_tciw_alt_lat)) * (gfdl4_tcc_alt_lat)


#---GFDL-AM4-AMIP Global Profile---#

gfdl4_tcc_alt_g = b['cf'][:] # 0-1
gfdl4_tclw_alt_g = b['lw'][:] #kg/kg
gfdl4_tciw_alt_g = b['iw'][:] #kg/kg
gfdl4_temp_alt_g = b['temp'][:] #K
gfdl4_plevel_alt_g = b['pressure'][:] #hPa

gfdl4_tcc_temp_g = b['cf_t'][:] # 0-1
gfdl4_tclw_temp_g = b['lw_t'][:] #kg/kg
gfdl4_tciw_temp_g = b['iw_t'][:] #kg/kg

gfdl4_tclw_frac_temp_g = b['lw_frac_t'][:]
gfdl4_tciw_frac_temp_g = b['iw_frac_t'][:]

#---GFDL-AM4-AMIP Phase Profile Fractions---#

gfdl4_tclw_frac_alt_g = b['lw_frac'][:]
gfdl4_tciw_frac_alt_g = b['iw_frac'][:]
gfdl4_tclw_frac_alt_so = b['lw_frac_so'][:]
gfdl4_tciw_frac_alt_so = b['iw_frac_so'][:]


#---GFDL-AM4-AMIP Southern Ocean Profile---#

gfdl4_tcc_alt_so = b['cf_so'][:23] # 0-1
gfdl4_tclw_alt_so = b['lw_so'][:19] #kg/kg
gfdl4_tciw_alt_so = b['iw_so'][:23] #kg/kg

gfdl4_tcc_temp_so = b['cf_t_so'][:] # 0-1
gfdl4_tclw_temp_so = b['lw_t_so'][:] #kg/kg
gfdl4_tciw_temp_so = b['iw_t_so'][:] #kg/kg

gfdl4_tclw_frac_temp_so = b['lw_frac_t_so'][:]
gfdl4_tciw_frac_temp_so = b['lw_frac_t_so'][:]


############################################################################### MRI_ESM2-AMIP Data

#---MRI_ESM2-AMIP Global Latitude Data---#

mri_tcc_lat_g = d['tcc'][:] # 0-1
mri_tclw_lat_g = d['tclw'][:] #kgm^-2
mri_tciw_lat_g = d['tciw'][:] #kgm^-2


#---MRI_ESM2-AMIP Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
mri_tclw_frac_lat_g = d['tclw_frac'][:] 
mri_tciw_frac_lat_g = d['tciw_frac'][:]



#---MRI_ESM2-AMIP Southern Ocean Latitude Data---#

mri_tcc_lat_so = mri_tcc_lat_g[mri_tcc_lat_g[:,0]>=-70]
mri_tcc_lat_so = mri_tcc_lat_so[mri_tcc_lat_so[:,0]<=-50] # 0-1

mri_tclw_lat_so = mri_tclw_lat_g[mri_tclw_lat_g[:,0]>=-70]
mri_tclw_lat_so = mri_tclw_lat_so[mri_tclw_lat_so[:,0]<=-50] #kgm^-2

mri_tciw_lat_so = mri_tciw_lat_g[mri_tciw_lat_g[:,0]>=-70]
mri_tciw_lat_so = mri_tciw_lat_so[mri_tciw_lat_so[:,0]<=-50] #kgm^-2


#---MRI_ESM2-AMIP lat-alt contour data---#

mri_tcc_alt_lat = d['cf_alt_lat'][:] #kg/kg
mri_tclw_alt_lat = d['lw_alt_lat'][:] #kg/kg
mri_tciw_alt_lat = d['iw_alt_lat'][:] #kg/kg
mri_temp_alt_lat = d['temp_alt_lat'][:] #kg/kg
mri_lat = d['lat'][:]
mri_alt = d['alt'][:]


mri_tclw_frac_alt_lat = (mri_tclw_alt_lat / (mri_tclw_alt_lat + mri_tciw_alt_lat)) * (mri_tcc_alt_lat / 100)
mri_tciw_frac_alt_lat = (mri_tciw_alt_lat / (mri_tclw_alt_lat + mri_tciw_alt_lat)) * (mri_tcc_alt_lat / 100)


#---MRI_ESM2-AMIP Global Profile---#

mri_tcc_alt_g = d['cf'][:] # 0-1
mri_tclw_alt_g = d['lw'][:] #kg/kg
mri_tciw_alt_g = d['iw'][:] #kg/kg
mri_temp_alt_g = d['temp'][:] #K
mri_plevel_alt_g = d['pressure'][:] #hPa

mri_tcc_temp_g = d['cf_t'][:] # 0-1
mri_tclw_temp_g = d['lw_t'][:] #kg/kg
mri_tciw_temp_g = d['iw_t'][:] #kg/kg

mri_tclw_frac_temp_g = d['lw_frac_t'][:]
mri_tciw_frac_temp_g = d['iw_frac_t'][:]


#---MRI_ESM2-AMIP Phase Profile Fractions---#

mri_tclw_frac_alt_g = d['lw_frac'][:]
mri_tciw_frac_alt_g = d['iw_frac'][:]
mri_tclw_frac_alt_so = d['lw_frac_so'][:]
mri_tciw_frac_alt_so = d['iw_frac_so'][:]


#---MRI_ESM2-AMIP Southern Ocean Profile---#

mri_tcc_alt_so = d['cf_so'][:] # 0-1
mri_tclw_alt_so = d['lw_so'][:] #kg/kg
mri_tciw_alt_so = d['iw_so'][:] #kg/kg

mri_tcc_temp_so = d['cf_t_so'][:] # 0-1
mri_tclw_temp_so = d['lw_t_so'][:] #kg/kg
mri_tciw_temp_so = d['iw_t_so'][:] #kg/kg

mri_tclw_frac_temp_so = d['lw_frac_t_so'][:]
mri_tciw_frac_temp_so = d['lw_frac_t_so'][:]



############################################################################### CESM2-CAM6-AMIP Data

#---CESM2-CAM6-AMIP Global Latitude Data---#

cam6_tcc_lat_g = e['tcc'][:] # 0-1
cam6_tclw_lat_g = e['tclw'][:] #kgm^-2
cam6_tciw_lat_g = e['tciw'][:] #kgm^-2


#---CESM2-CAM6-AMIP Global Latitude Phase Fractions---#

# 0-1 derived from fraction of mean liquid and ice water content at specific latitude * total cloud fraction at the latitude
cam6_tclw_frac_lat_g = e['tclw_frac'][:] 
cam6_tciw_frac_lat_g = e['tciw_frac'][:]


#---CESM2-CAM6-AMIP Southern Ocean Latitude Data---#

cam6_tcc_lat_so = cam6_tcc_lat_g[cam6_tcc_lat_g[:,0]>=-70]
cam6_tcc_lat_so = cam6_tcc_lat_so[cam6_tcc_lat_so[:,0]<=-50] # 0-1

cam6_tclw_lat_so = cam6_tclw_lat_g[cam6_tclw_lat_g[:,0]>=-70]
cam6_tclw_lat_so = cam6_tclw_lat_so[cam6_tclw_lat_so[:,0]<=-50] #kgm^-2

cam6_tciw_lat_so = cam6_tciw_lat_g[cam6_tciw_lat_g[:,0]>=-70]
cam6_tciw_lat_so = cam6_tciw_lat_so[cam6_tciw_lat_so[:,0]<=-50] #kgm^-2


#---CESM2-CAM6-AMIP lat-alt contour data---#

cam6_tcc_alt_lat = e['cf_alt_lat'][:] 
cam6_tclw_alt_lat = e['lw_alt_lat'][:] #kg/kg
cam6_tciw_alt_lat = e['iw_alt_lat'][:] #kg/kg
cam6_temp_alt_lat = e['temp_alt_lat'][:] #kg/kg
cam6_lat = e['lat'][:]
cam6_lat = np.hstack(cam6_lat)
cam6_alt = e['alt'][:]
cam6_alt = np.hstack(cam6_alt)

cam6_tclw_frac_alt_lat = (cam6_tclw_alt_lat / (cam6_tclw_alt_lat + cam6_tciw_alt_lat)) * (cam6_tcc_alt_lat / 100)
cam6_tciw_frac_alt_lat = (cam6_tciw_alt_lat / (cam6_tclw_alt_lat + cam6_tciw_alt_lat)) * (cam6_tcc_alt_lat / 100)


#---CESM2-CAM6-AMIP Global Profile---#

cam6_tcc_alt_g = e['cf'][:] # 0-1
cam6_tclw_alt_g = e['lw'][:] #kg/kg
cam6_tciw_alt_g = e['iw'][:] #kg/kg
cam6_temp_alt_g = e['temp'][:] #K
cam6_plevel_alt_g = e['pressure'][:] #hPa

cam6_tcc_temp_g = e['cf_t'][:] # 0-1
cam6_tclw_temp_g = e['lw_t'][:] #kg/kg
cam6_tciw_temp_g = e['iw_t'][:] #kg/kg

cam6_tclw_frac_temp_g = e['lw_frac_t'][:]
cam6_tciw_frac_temp_g = e['iw_frac_t'][:]


#---CESM2-CAM6-AMIP Phase Profile Fractions---#

cam6_tclw_frac_alt_g = e['lw_frac'][:]
cam6_tciw_frac_alt_g = e['iw_frac'][:]
cam6_tclw_frac_alt_so = e['lw_frac_so'][:]
cam6_tciw_frac_alt_so = e['iw_frac_so'][:]


#---CESM2-CAM6-AMIP Southern Ocean Profile---#

cam6_tcc_alt_so = e['cf_so'][:] # 0-1
cam6_tclw_alt_so = e['lw_so'][:] #kg/kg
cam6_tciw_alt_so = e['iw_so'][:] #kg/kg

cam6_tcc_temp_so = e['cf_t_so'][:] # 0-1
cam6_tclw_temp_so = e['lw_t_so'][:] #kg/kg
cam6_tciw_temp_so = e['iw_t_so'][:] #kg/kg

cam6_tclw_frac_temp_so = e['lw_frac_t_so'][:]
cam6_tciw_frac_temp_so = e['lw_frac_t_so'][:]


############################################################################### CALIPSO Data

#---CALIPSO Global Latitude Data---#

calipso_tcc_lat_g = f['tcc'][:] # 0-1
calipso_tclw_frac_lat_g = f['tclw_frac'][:] # 0-1
calipso_tciw_frac_lat_g = f['tciw_frac'][:] # 0-1

#---CALIPSO Southern Ocean Latitude Data---#

calipso_tcc_lat_so = calipso_tcc_lat_g[calipso_tcc_lat_g[:,0]>=-70]
calipso_tcc_lat_so = calipso_tcc_lat_so[calipso_tcc_lat_so[:,0]<=-50] # 0-1

calipso_tclw_frac_lat_so = calipso_tclw_frac_lat_g[calipso_tclw_frac_lat_g[:,0]>=-70]
calipso_tclw_frac_lat_so = calipso_tclw_frac_lat_so[calipso_tclw_frac_lat_so[:,0]<=-50] # 0-1

calipso_tciw_frac_lat_so = calipso_tciw_frac_lat_g[calipso_tciw_frac_lat_g[:,0]>=-70]
calipso_tciw_frac_lat_so = calipso_tciw_frac_lat_so[calipso_tciw_frac_lat_so[:,0]<=-50] # 0-1


#---CALIPSO Global Profile---#

calipso_tcc_alt_g = f['cf'][:] # 0-1
calipso_tclw_frac_alt_g = f['lw_frac'][:] # 0-1
calipso_tciw_frac_alt_g = f['iw_frac'][:] # 0-1


#---CALIPSO lat-alt contour data---#

calipso_tcc_alt_lat = f['cf_alt_lat'][:]
calipso_tclw_frac_alt_lat = f['liq_frac_alt_lat'][:]
calipso_tciw_frac_alt_lat = f['ice_frac_alt_lat'][:]
calipso_lat = f['lat'][:]
calipso_lat = np.hstack(calipso_lat)
calipso_alt = f['alt'][:]
calipso_alt = np.hstack(calipso_alt)
calipso_temp = f['alt_t'][:]
calipso_temp = np.hstack(calipso_temp)

#---CALIPSO lat-temp contour data---#

calipso_cf_t_lat = f['cf_t_lat'][:]
calipso_lw_t_lat = f['lw_t_lat'][:]
calipso_iw_t_lat = f['iw_t_lat'][:]

#---CALIPSO temp profile data---#

calipso_tcc_temp_g = f['cf_t'][:] # 0-1
calipso_tclw_frac_temp_g = f['lw_t_frac'][:] # 0-1
calipso_tciw_frac_temp_g = f['iw_t_frac'][:] # 0-1


#---CALIPSO Southern Ocean Profile---#


calipso_tcc_alt_so = f['cf_so'][:] # 0-1
calipso_tclw_frac_alt_so = f['lw_frac_so'][:] # 0-1
calipso_tciw_frac_alt_so = f['iw_frac_so'][:] # 0-1

calipso_tcc_temp_so = f['cf_t_so'][:] # 0-1
calipso_tclw_frac_temp_so = f['lw_t_frac_so'][:] # 0-1
calipso_tciw_frac_temp_so = f['iw_t_frac_so'][:] # 0-1

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


"""

os.chdir('E:/University/University/MSc/Models/Images/Meeting 1.7')

#os.chdir('C:/Users/toha006/University/University/MSc/Models/Images/Meeting 1.7')
"""

############################################################################### Temperature Profiles

#---Plot Satellite Liquid Cloud Fraction with Temperature---#

"""

fig, ax = plt.subplots()


ax.plot(cccm_tclw_frac_temp_g[5:44,0],cccm_tclw_frac_temp_g[5:44,1], '-r', label='CCCM - Global')
ax.plot(cccm_tclw_frac_temp_so[2:40,0],cccm_tclw_frac_temp_so[2:40,1], '--r', label='CCCM - Southern Ocean')

ax.plot(calipso_tclw_frac_temp_g[20:34,0],calipso_tclw_frac_temp_g[20:34,1], '-b', label='CALIPSO - Global')
ax.plot(calipso_tclw_frac_temp_so[20:34,0],calipso_tclw_frac_temp_so[20:34,1], '--b', label='CALIPSO - Southern Ocean')
plt.axvline(x=273, label = '273K', color = 'green')
ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Liquid Cloud Fraction')

plt.title('07.2006 to 04.2011 Satellite Liquid Cloud Fraction vs Temperature')
plt.grid(True)
plt.savefig("07.2006_04.2011_satellite_liquid_T.svg", format="svg", bbox_inches='tight')

plt.show()
"""


#---Plot Liquid Cloud Fraction with Temperature---#


"""
fig, ax = plt.subplots()

ax.plot(cccm_tclw_frac_temp_g[5:44,0],cccm_tclw_frac_temp_g[5:44,1], '-r', label='Global - CCCM')
ax.plot(calipso_tclw_frac_temp_g[20:34,0],calipso_tclw_frac_temp_g[20:34,1], '-b', label='Global - CALIPSO')
ax.plot(ecmwf_tclw_frac_temp_g[18:,0],ecmwf_tclw_frac_temp_g[18:,1], '-k', label='Global - ECMWF-ERA5')
ax.plot(gfdl4_tclw_frac_temp_g[:18,0],gfdl4_tclw_frac_temp_g[:18,1], '-g', label='Global - CMIP6-GFDL-AM4-AMIP')
ax.plot(mri_tclw_frac_temp_g[:26,0],mri_tclw_frac_temp_g[:26,1], '-m', label='Global - CMIP6-MRI_ESM2-AMIP')
ax.plot(cam6_tclw_frac_temp_g[17:,0],cam6_tclw_frac_temp_g[17:,1], '-c', label='Global - CMIP6-CESM2-CAM6-AMIP')

ax.plot(cccm_tclw_frac_temp_so[2:40,0],cccm_tclw_frac_temp_so[2:40,1], '--r', label='Southern Ocean - CCCM')
ax.plot(calipso_tclw_frac_temp_so[20:34,0],calipso_tclw_frac_temp_so[20:34,1], '--b', label='Southern Ocean - CALIPSO')
ax.plot(ecmwf_tclw_frac_temp_so[20:,0],ecmwf_tclw_frac_temp_so[20:,1], '--k', label='Southern Ocean - ECMWF-ERA5')
ax.plot(gfdl4_tclw_frac_temp_so[:18,0],gfdl4_tclw_frac_temp_so[:18,1], '--g', label='Southern Ocean - CMIP6-GFDL-AM4-AMIP')
ax.plot(mri_tclw_frac_temp_so[:26,0],mri_tclw_frac_temp_so[:26,1], '--m', label='Southern Ocean - CMIP6-MRI_ESM2-AMIP')
ax.plot(cam6_tclw_frac_temp_so[17:,0],cam6_tclw_frac_temp_so[17:,1], '--c', label='Southern Ocean - CMIP6-CESM2-CAM6-AMIP')


plt.axvline(x=273, label = '273K', color = 'red')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Liquid Cloud Fraction')

plt.title('07.2006 to 04.2011 Liquid Cloud Fraction vs Temperature')

plt.grid(True)
plt.savefig("07.2006_04.2011_liquid_T.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Liquid Cloud Fraction with Temperature---#


"""
fig, ax = plt.subplots()

ax.plot(cccm_tclw_temp_g[20:36,0],cccm_tclw_temp_g[20:36,1], '-b', label='Global - CALIPSO')
ax.plot(ecmwf_tclw_temp_g[18:,0],ecmwf_tclw_temp_g[18:,1], '-k', label='Global - ECMWF-ERA5')
ax.plot(gfdl4_tclw_temp_g[:18,0],gfdl4_tclw_temp_g[:18,1], '-g', label='Global - CMIP6-GFDL-AM4')
ax.plot(mri_tclw_temp_g[:26,0],mri_tclw_temp_g[:26,1], '-m', label='Global - MRI_ESM2-AMIP')
ax.plot(cam6_tclw_temp_g[17:,0],cam6_tclw_temp_g[17:,1], '-c', label='Global - CESM2-CAM6-AMIP')

ax.plot(cccm_tclw_temp_so[20:36,0],cccm_tclw_temp_so[20:36,1], '--b', label='Southern Ocean - CALIPSO')
ax.plot(ecmwf_tclw_temp_so[20:,0],ecmwf_tclw_temp_so[20:,1], '--k', label='Southern Ocean - ECMWF-ERA5')
ax.plot(gfdl4_tclw_temp_so[:18,0],gfdl4_tclw_temp_so[:18,1], '--g', label='Southern Ocean - CMIP6-GFDL-AM4')
ax.plot(mri_tclw_temp_so[:26,0],mri_tclw_temp_so[:26,1], '--m', label='Southern Ocean - MRI_ESM2-AMIP')
ax.plot(cam6_tclw_temp_so[17:,0],cam6_tclw_temp_so[17:,1], '--c', label='Southern Ocean - CESM2-CAM6-AMIP')


plt.axvline(x=273, label = '273K', color = 'red')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Liquid Cloud Fraction')

plt.title('07.2006 to 04.2011 Liquid Cloud Fraction vs Temperature')

plt.grid(True)
plt.savefig("07.2006_04.2011_liquid_T.svg", format="svg", bbox_inches='tight')
plt.show()
"""


############################################################################### Global Latitude Plots

#---Plot Global Cloud Fraction with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tcc_lat_g[:,0],cccm_tcc_lat_g[:,1], '-r', label='CCCM')
ax.plot(calipso_tcc_lat_g[:,0],calipso_tcc_lat_g[:,1], '-b', label='CAPLISO')
ax.plot(ceres_tcc_lat_g[:,0],ceres_tcc_lat_g[:,1], '--r', label='CERES')

ax.plot(ecmwf_tcc_lat_g[:,0],ecmwf_tcc_lat_g[:,1], '-k', label='ECMWF-ERA5')
ax.plot(gfdl4_tcc_lat_g[:,0],gfdl4_tcc_lat_g[:,1], '-g', label='CMIP6-GFDL-AM4-AMIP')
ax.plot(mri_tcc_lat_g[:,0],mri_tcc_lat_g[:,1], '-m', label='CMIP6-MRI_ESM2-AMIP')
ax.plot(cam6_tcc_lat_g[:,0],cam6_tcc_lat_g[:,1], '-c', label='CMIP6-CESM2-CAM6-AMIP')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 Global Cloud Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tcc_lat_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""


#---Plot Global Comparison Liquid Cloud Fractions with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tclw_frac_lat_g[:,0],cccm_tclw_frac_lat_g[:,1], ':r', label='CCCM')
ax.plot(calipso_tciw_frac_lat_g[:,0],calipso_tciw_frac_lat_g[:,1], ':b', label='CALIPSO')
ax.plot(ecmwf_tclw_frac_lat_g[:,0],ecmwf_tclw_frac_lat_g[:,1], ':k', label='ECMWF-ERA5')
ax.plot(gfdl4_tclw_frac_lat_g[:,0],gfdl4_tclw_frac_lat_g[:,1], ':g', label='CMIP6-GFDL-AM4-AMIP')
ax.plot(mri_tclw_frac_lat_g[:,0],mri_tclw_frac_lat_g[:,1], ':m', label='CMIP6-MRI_ESM2-AMIP')
ax.plot(cam6_tclw_frac_lat_g[:,0],cam6_tclw_frac_lat_g[:,1], ':c', label='CMIP6-CESM2-CAM6-AMIP')

ax.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Liquid Water Fraction')

plt.title('07.2006 to 04.2011 Global Liquid Fractions of Clouds vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot Global Comparison Ice Cloud Fractions with Latitude---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(cccm_tciw_frac_lat_g[:,0],cccm_tciw_frac_lat_g[:,1], '--r', label='CCCM')
ax.plot(calipso_tclw_frac_lat_g[:,0],calipso_tclw_frac_lat_g[:,1], '--b', label='CALIPSO')
ax.plot(ecmwf_tciw_frac_lat_g[:,0],ecmwf_tciw_frac_lat_g[:,1], '--k', label='ECMWF-ERA5')
ax.plot(gfdl4_tciw_frac_lat_g[:,0],gfdl4_tciw_frac_lat_g[:,1], '--g', label='CMIP6-GFDL-AM4')
ax.plot(mri_tciw_frac_lat_g[:,0],mri_tciw_frac_lat_g[:,1], '--m', label='MRI_ESM2-AMIP')
ax.plot(cam6_tciw_frac_lat_g[:,0],cam6_tciw_frac_lat_g[:,1], '--c', label='CESM2-CAM6-AMIP')


ax.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax.set_xlabel('Latitude')
ax.set_ylabel('Cloud Ice Water Fraction')

plt.title('07.2006 to 04.2011 Global Ice Fractions of Clouds vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tciw_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CALIPSO Cloud Fraction and Phase Fraction with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(calipso_tcc_lat_g[:,0],calipso_tcc_lat_g[:,1], '-b', label='Cloud Fraction')
ax1.plot(calipso_tciw_frac_lat_g[:,0],calipso_tciw_frac_lat_g[:,1], ':b', label='Liquid Fraction')
ax1.plot(calipso_tclw_frac_lat_g[:,0],calipso_tclw_frac_lat_g[:,1], '--b', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(1.4, 1.0));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CALIPSO Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_CALIPSO_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot ECMWF Cloud Fraction and Phase Content with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(ecmwf_tcc_lat_g[:,0],ecmwf_tcc_lat_g[:,1], '-k', label='Cloud Fraction')
ax1.plot(ecmwf_tclw_frac_lat_g[:,0],ecmwf_tclw_frac_lat_g[:,1], ':k', label='Liquid Fraction')
ax1.plot(ecmwf_tciw_frac_lat_g[:,0],ecmwf_tciw_frac_lat_g[:,1], '--k', label='Ice Fraction')

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

ax1.plot(gfdl4_tcc_lat_g[:,0],gfdl4_tcc_lat_g[:,1], '-g', label='Cloud Fraction')
ax1.plot(gfdl4_tclw_frac_lat_g[:,0],gfdl4_tclw_frac_lat_g[:,1], ':g', label='Liquid Fraction')
ax1.plot(gfdl4_tciw_frac_lat_g[:,0],gfdl4_tciw_frac_lat_g[:,1], '--g', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CMIP6-GFDL-AM4 Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_gfdl4_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""


#---Plot MRI_ESM2-AMIP Cloud Fraction and Phase Fraction with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(mri_tcc_lat_g[:,0],mri_tcc_lat_g[:,1], '-m', label='Cloud Fraction')
ax1.plot(mri_tclw_frac_lat_g[:,0],mri_tclw_frac_lat_g[:,1], ':m', label='Liquid Fraction')
ax1.plot(mri_tciw_frac_lat_g[:,0],mri_tciw_frac_lat_g[:,1], '--m', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 MRI_ESM2-AMIP Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_MRI_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CAM6 Cloud Fraction and Phase Fraction with Latitude---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(cam6_tcc_lat_g[:,0],cam6_tcc_lat_g[:,1], '-c', label='Cloud Fraction')
ax1.plot(cam6_tclw_frac_lat_g[:,0],cam6_tclw_frac_lat_g[:,1], ':c', label='Liquid Fraction')
ax1.plot(cam6_tciw_frac_lat_g[:,0],cam6_tciw_frac_lat_g[:,1], '--c', label='Ice Fraction')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_xlabel('Latitude')
ax1.set_ylabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CMIP6-CESM2.1-CAM6 Global Cloud Fraction and Phase Fraction vs Latitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_CAM6_lat_frac_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""


############################################################################### Global Altitude Plots

#---Plot Global Cloud Fraction Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(calipso_tcc_alt_g[:,1],calipso_tcc_alt_g[:,0], '-b', label='CALIPSO')

ax.plot(ecmwf_tcc_alt_g[:,1],ecmwf_tcc_alt_g[:,0], '-k', label='ECMWF-ERA5')
ax.plot(gfdl4_tcc_alt_g[:,1],gfdl4_tcc_alt_g[:,0], '-g', label='GFDL-AM4-AMIP')
ax.plot(mri_tcc_alt_g[:,1],mri_tcc_alt_g[:,0], '-m', label='MRI-ESM2-AMIP')
ax.plot(cam6_tcc_alt_g[:,1],cam6_tcc_alt_g[:,0], '-c', label='CESM2.1-CAM6-AMIP')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 Global Cloud Fraction vs Altitude')

plt.grid(True)

plt.savefig("07.2006_04.2011_tcc_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Cloud Liquid Water Fraction Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(calipso_tclw_frac_alt_g[:21,1],calipso_tclw_frac_alt_g[:21,0], '-b', label='CALIPSO')
ax.plot(ecmwf_tclw_frac_alt_g[15:,1],ecmwf_tclw_frac_alt_g[15:,0], '-k', label='ECMWF-ERA5')
ax.plot(gfdl4_tclw_frac_alt_g[:20,1],gfdl4_tclw_frac_alt_g[:20,0], '-g', label='GFDL-AM4-AMIP')
ax.plot(mri_tclw_frac_alt_g[:30,1],mri_tclw_frac_alt_g[:30,0], '-m', label='MRI-ESM2-AMIP')
ax.plot(cam6_tclw_frac_alt_g[15:,1],cam6_tclw_frac_alt_g[15:,0], '-c', label='CESM2.1-CAM6-AMIP')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Cloud Liquid Water Fraction')

plt.title('07.2006 to 04.2011 Cloud Liquid Water Fraction vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_frac_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Global Cloud Ice Water Fraction Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(calipso_tciw_frac_alt_g[:,1],calipso_tciw_frac_alt_g[:,0], '--b', label='CALIPSO')
ax.plot(ecmwf_tciw_frac_alt_g[10:,1],ecmwf_tciw_frac_alt_g[10:,0], '--k', label='ECMWF-ERA5')
ax.plot(gfdl4_tciw_frac_alt_g[:26,1],gfdl4_tciw_frac_alt_g[:26,0], '--g', label='GFDL-AM4-AMIP')
ax.plot(mri_tciw_frac_alt_g[:50,1],mri_tciw_frac_alt_g[:50,0], '--m', label='MRI-ESM2-AMIP')
ax.plot(cam6_tciw_frac_alt_g[7:,1],cam6_tciw_frac_alt_g[7:,0], '--c', label='CESM2.1-CAM6-AMIP')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Cloud Ice Water Fraction')

plt.title('07.2006 to 04.2011 Cloud Ice Water Fraction vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tciw_frac_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""


#---Plot Global Specific Liquid Water Content with Temperature Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(ecmwf_tclw_temp_g[:,1],ecmwf_tclw_temp_g[:,0], ':k', label='ECMWF-ERA5')
ax.plot(gfdl4_tclw_temp_g[:,1],gfdl4_tclw_temp_g[:,0], ':g', label='GFDL-AM4-AMIP')
ax.plot(mri_tclw_temp_g[:,1],mri_tclw_temp_g[:,0], ':m', label='MRI-ESM2-AMIP')
ax.plot(cam6_tclw_temp_g[:,1],cam6_tclw_temp_g[:,0], ':c', label='CESM2.1-CAM6-AMIP')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

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

ax.plot(ecmwf_tciw_temp_g[:,1]*10000,ecmwf_tciw_temp_g[:,0], '-k', label='ECMWF-ERA5')
ax.plot(gfdl3_tciw_temp_g[:,1]*10000,gfdl3_tciw_temp_g[:,0], '-b', label='CMIP5-GFDL-CM3')
ax.plot(gfdl4_tciw_temp_g[:,1]*10000,gfdl4_tciw_temp_g[:,0], '-g', label='CMIP6-GFDL-AM4')
ax.plot(cam_tciw_temp_g[:,1]*10000,cam_tciw_temp_g[:,0], '-m', label='CMIP5-CESM1-CAM5-RPC4.5')
ax.plot(cam6_tciw_temp_g[:,1]*10000,cam6_tciw_temp_g[:,0], '-c', label='CMIP6-CESM2.1-CAM6-RPC4.5')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

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

ax.plot(ecmwf_tclw_temp_g[:,1]*10000,ecmwf_tclw_temp_g[:,0], '-k', label='ECMWF-ERA5')
ax.plot(gfdl3_tclw_temp_g[:,1]*10000,gfdl3_tclw_temp_g[:,0], '-b', label='CMIP5-GFDL-CM3')
ax.plot(gfdl4_tclw_temp_g[:,1]*10000,gfdl4_tclw_temp_g[:,0], '-g', label='CMIP6-GFDL-AM4')
ax.plot(cam_tclw_temp_g[:,1]*10000,cam_tclw_temp_g[:,0], '-m', label='CMIP5-CESM1-CAM5-RPC4.5')
ax.plot(cam6_tclw_temp_g[:,1]*10000,cam6_tclw_temp_g[:,0], '-c', label='CMIP6-CESM2.1-CAM6-RPC4.5')

ax.plot(ecmwf_tciw_temp_g[:,1]*10000,ecmwf_tciw_temp_g[:,0], '-k', label='ECMWF-ERA5')
ax.plot(gfdl3_tciw_temp_g[:,1]*10000,gfdl3_tciw_temp_g[:,0], '-b', label='CMIP5-GFDL-CM3')
ax.plot(gfdl4_tciw_temp_g[:,1]*10000,gfdl4_tciw_temp_g[:,0], '-g', label='CMIP6-GFDL-AM4')
ax.plot(cam_tciw_temp_g[:,1]*10000,cam_tciw_temp_g[:,0], '-m', label='CMIP5-CESM1-CAM5-RPC4.5')
ax.plot(cam6_tciw_temp_g[:,1]*10000,cam6_tciw_temp_g[:,0], '-c', label='CMIP6-CESM2.1-CAM6-RPC4.5')


ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_ylabel('Temperature (K)')
ax.set_xlabel('Specific Content $(kg/kg) x 10^{-4}$')

plt.title('07.2006 to 04.2011 Global Specific Liquid and Ice Water Content vs Temperature Profile')
plt.gca().invert_yaxis()
plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_tciw_T_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot CALIPSO Cloud and Phase Fraction Profile---#

"""
plt.figure()
fig, ax1 = plt.subplots()

ax1.plot(calipso_tcc_alt_g[:,1],calipso_tcc_alt_g[:,0], '-b', label='Cloud Fraction')
ax1.plot(calipso_tclw_frac_alt_g[:,1],calipso_tclw_frac_alt_g[:,0], ':b', label='Liquid Content')
ax1.plot(calipso_tciw_frac_alt_g[:,1],calipso_tciw_frac_alt_g[:,0], '--b', label='Ice Content')

ax1.legend(loc='upper center', bbox_to_anchor=(0.3, -0.15));

ax1.set_ylabel('Altitude')
ax1.set_xlabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 CALIPSO Global Cloud Fraction and Phase Content vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_CALIPSO_alt_g.svg", format="svg", bbox_inches='tight')
plt.show()
"""

############################################################################### Southern Ocean Altitude Plots

#---Plot Southern Ocean Cloud Fraction Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(calipso_tcc_alt_so[:,1],calipso_tcc_alt_so[:,0], '-r', label='CALIPSO')

ax.plot(ecmwf_tcc_alt_so[:,1],ecmwf_tcc_alt_so[:,0], '-k', label='ECMWF-ERA5')
ax.plot(gfdl3_tcc_alt_so[:,1],gfdl3_tcc_alt_so[:,0], '-b', label='CMIP5-GFDL-CM3')
ax.plot(gfdl4_tcc_alt_so[:,1],gfdl4_tcc_alt_so[:,0], '-g', label='CMIP6-GFDL-AM4')
ax.plot(cam_tcc_alt_so[:,1],cam_tcc_alt_so[:,0], '-m', label='CMIP5-CESM1-CAM5-RPC4.5')
ax.plot(cam6_tcc_alt_so[:,1],cam6_tcc_alt_so[:,0], '-c', label='CMIP6-CESM2.1-CAM6-RPC4.5')

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Cloud Fraction')

plt.title('07.2006 to 04.2011 Southern Ocean Cloud Fraction vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tcc_alt_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Cloud Liquid Water Fraction Altitude Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(calipso_tclw_frac_alt_so[:21,1],calipso_tclw_frac_alt_so[:21,0], '-r', label='CALIPSO')
ax.plot(ecmwf_tclw_frac_alt_so[17:36,1],ecmwf_tclw_frac_alt_so[17:36,0], '-k', label='ECMWF-ERA5')
ax.plot(gfdl3_tclw_frac_alt_so[:20,1],gfdl3_tclw_frac_alt_so[:20,0], '-b', label='CMIP5-GFDL-CM3')
ax.plot(gfdl4_tclw_frac_alt_so[:20,1],gfdl4_tclw_frac_alt_so[:20,0], '-g', label='CMIP6-GFDL-AM4')
ax.plot(cam_tclw_frac_alt_so[16:,1],cam_tclw_frac_alt_so[16:,0], '-m', label='CMIP5-CESM1-CAM5-RPC4.5')
ax.plot(cam6_tclw_frac_alt_so[16:,1],cam6_tclw_frac_alt_so[16:,0], '-c', label='CMIP6-CESM2.1-CAM6-RPC4.5')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Cloud Liquid Water Fraction')

plt.title('07.2006 to 04.2011 Cloud Liquid Water Fraction over the Southern Ocean vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tclw_frac_alt_so.svg", format="svg", bbox_inches='tight')

plt.show()
"""

#---Plot Southern Ocean Cloud Ice Water Content Fraction Profile---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.plot(calipso_tciw_frac_alt_so[:,1],calipso_tciw_frac_alt_so[:,0], '-r', label='CALIPSO')
ax.plot(ecmwf_tciw_frac_alt_so[:,1],ecmwf_tciw_frac_alt_so[:,0], '-k', label='ECMWF-ERA5')
ax.plot(gfdl3_tciw_frac_alt_so[:,1],gfdl3_tciw_frac_alt_so[:,0], '-b', label='CMIP5-GFDL-CM3')
ax.plot(gfdl4_tciw_frac_alt_so[:,1],gfdl4_tciw_frac_alt_so[:,0], '-g', label='CMIP6-GFDL-AM4')
ax.plot(cam_tciw_frac_alt_so[:,1],cam_tciw_frac_alt_so[:,0], '-m', label='CMIP5-CESM1-CAM5-RPC4.5')
ax.plot(cam6_tciw_frac_alt_so[:,1],cam6_tciw_frac_alt_so[:,0], '-c', label='CMIP6-CESM2.1-CAM6-RPC4.5')

ax.legend(loc='upper center', bbox_to_anchor=(1.3, 1.0));

ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Cloud Ice Water Fraction')

plt.title('07.2006 to 04.2011 Cloud Ice Water Fraction over the Southern Ocean vs Altitude')

plt.grid(True)
plt.savefig("07.2006_04.2011_tciw_frac_alt_g.svg", format="svg", bbox_inches='tight')

plt.show()
"""
############################################################################### Contour plots

#---Plot CALIPSO Cloud Liquid Water Fraction---#

"""
plt.figure()
fig, ax = plt.subplots()

cs = ax.contourf(calipso_lat, calipso_alt[:16], calipso_tclw_frac_alt_lat[:16])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')

cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Liquid Cloud Fraction')

plt.title('07.2006 to 04.2011 CALIPSO Cloud Liquid Water Fraction')

plt.savefig("07.2006_04.2011_contour_calipso_tclw.svg", format="svg", bbox_inches='tight')
plt.show()

"""

#---Plot CALIPSO Cloud Ice Water Fraction---#

"""
plt.figure()
fig, ax = plt.subplots()

cs = ax.contourf(calipso_lat, calipso_alt, calipso_tciw_frac_alt_lat)

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Ice Cloud Fraction')

plt.title('07.2006 to 04.2011 CALIPSO Cloud Ice Water Fraction')

plt.savefig("07.2006_04.2011_contour_calipso_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot ECMWF Cloud Liquid Water Fraction---#

"""
plt.figure()
fig, ax = plt.subplots()

cs = ax.contourf(ecmwf_lat, ecmwf_alt[0:18], ecmwf_tclw_frac_alt_lat[0:18])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Liquid Cloud Fraction')

plt.title('07.2006 to 04.2011 ECMWF-ERA5 Cloud Liquid Water Fraction')

plt.savefig("07.2006_04.2011_contour_ecmwf_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot ECMWF Cloud Ice Water Fraction---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.contourf(ecmwf_lat, ecmwf_alt[:26], ecmwf_tciw_frac_alt_lat[:26])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Ice Cloud Fraction')

plt.title('07.2006 to 04.2011 ECMWF-ERA5 Cloud Ice Water Fraction')

plt.savefig("07.2006_04.2011_contour_ecmwf_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""


#---Plot GFDL-AM4-AMIP Cloud Liquid Water Fraction---#

"""
plt.figure()
fig, ax = plt.subplots()

cs = ax.contourf(gfdl4_lat, gfdl4_alt[0:18], gfdl4_tclw_frac_alt_lat[0:18])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Liquid Cloud Fraction')

plt.title('07.2006 to 04.2011 GFDL-AM4-AMIP Cloud Liquid Water Fraction')

plt.savefig("07.2006_04.2011_contour_gfdl4_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot GFDL-AM4-AMIP Cloud Ice Water Fraction---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.contourf(gfdl4_lat, gfdl4_alt[:26], gfdl4_tciw_frac_alt_lat[:26])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Ice Cloud Fraction')

plt.title('07.2006 to 04.2011 GFDL-AM4-AMIP Cloud Ice Water Fraction')

plt.savefig("07.2006_04.2011_contour_gfdl4_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""
#---Plot MRI Cloud Liquid Water Fraction---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.contourf(mri_lat, mri_alt[0:25], mri_tclw_frac_alt_lat[0:25])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Liquid Cloud Fraction')

plt.title('07.2006 to 04.2011 MRI_ESM2-AMIP Cloud Liquid Water Fraction')

plt.savefig("07.2006_04.2011_contour_mri_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot MRI Cloud Ice Water Fraction---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.contourf(mri_lat, mri_alt[:43], mri_tciw_frac_alt_lat[:43])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Ice Cloud Fraction')

plt.title('07.2006 to 04.2011 MRI_ESM2-AMIP Cloud Ice Water Fraction')

plt.savefig("07.2006_04.2011_contour_mri_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CESM2-CAM6-AMIP Cloud Liquid Water Fraction---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.contourf(cam6_lat, cam6_alt[18:32], cam6_tclw_frac_alt_lat[18:32])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Liquid Cloud Fraction')

plt.title('07.2006 to 04.2011 CESM2-CAM6-AMIP Cloud Liquid Water Fraction')

plt.savefig("07.2006_04.2011_contour_CAM6_tclw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Plot CESM2-CAM6-AMIP Cloud Ice Water Fraction---#

"""
plt.figure()
fig, ax = plt.subplots()

ax.contourf(cam6_lat, cam6_alt[7:32], cam6_tciw_frac_alt_lat[7:32])

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);

ax.set_xlabel('Latitude')
ax.set_ylabel('Altitude (km)')
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Ice Cloud Fraction')

plt.title('07.2006 to 04.2011 CESM2-CAM6-AMIP Cloud Ice Water Fraction')

plt.savefig("07.2006_04.2011_contour_CAM6_tciw.svg", format="svg", bbox_inches='tight')
plt.show()
"""

#---Sub-Plot CALIPSO and GFDL4 Cloud Liquid Water Fraction---#

"""
plt.figure()

plt.subplot(1, 2, 1)
plt.contourf(calipso_lat, calipso_alt[:16], calipso_tclw_alt_lat[:16])
plt.xlabel('Latitude')
plt.ylabel('Altitude (km)')
plt.subplot(1, 2, 2)
plt.contourf(calipso_lat, calipso_alt[:16], calipso_tclw_alt_lat[:16])

plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4);



plt.title('07.2006 to 04.2011 CALIPSO Cloud Liquid Water Fraction')

#plt.savefig("07.2006_04.2011_contour_calipso_tclw.svg", format="svg", bbox_inches='tight')
plt.show()

"""